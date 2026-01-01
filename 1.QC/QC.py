import pandas as pd
from pyliftover import LiftOver
from tqdm import tqdm
import os
import glob
import numpy as np

# ================= 配置区域 =================
# 1. 输入数据的文件夹路径
INPUT_FOLDER = "/home/chenguanglei2011/article/0.Raw/ukb"

# 2. 指定输入数据的版本 ('hg38' 或 'hg19')
INPUT_BUILD = 'hg38'  # <--- 请根据你的数据来源修改这里

# 3. 定义 MHC 区域 (固定基于 hg19 坐标)
MHC_CHROM = '6'  
MHC_START = 28477797
MHC_END = 33448354

# 4. MAF 阈值
MAF_THRESHOLD = 0.01

# 5. 需要检查空值的关键列
CRITICAL_COLS = ['SNP', 'effect_allele', 'other_allele', 'eaf', 'beta', 'se', 'pval', 'chr', 'pos']
# ===========================================

def main():
    print(f"当前设置: 输入数据版本为 【{INPUT_BUILD}】")
    
    dirs = ["1.hg19", "2.hg38", "3.FUMA"]
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    print("正在初始化转换链...")
    try:
        lo_38_to_19 = LiftOver('hg38', 'hg19')
        lo_19_to_38 = LiftOver('hg19', 'hg38')
    except Exception as e:
        print(f"初始化 LiftOver 失败: {e}")
        return

    files = glob.glob(os.path.join(INPUT_FOLDER, "*.txt"))
    print(f"共发现 {len(files)} 个文件需要处理。")

    for file_path in files:
        file_name = os.path.basename(file_path)
        
        # ==================================================
        # 【新增】断点续传检测
        # ==================================================
        # 检查最终输出目录(2.hg38)是否已有该文件
        final_output_path = os.path.join("2.hg38", file_name)
        if os.path.exists(final_output_path):
            print(f">>> [跳过] 文件已存在 (断点续传): {file_name}")
            continue
        # ==================================================

        print(f"\n{'='*50}")
        print(f"正在处理文件: {file_name} (输入视作 {INPUT_BUILD})")
        
        # --- 读取数据 ---
        try:
            df = pd.read_csv(file_path, sep=r'\s+', low_memory=False)
        except Exception as e:
            print(f"读取失败，跳过: {e}")
            continue

        raw_len = len(df)
        
        # ==================================================
        # 深度数据清洗
        # ==================================================
        print(f">> 数据清洗前行数: {raw_len}")
        
        # 1. 强制转数值 (处理可能存在的字符型数字)
        numeric_cols = ['eaf', 'beta', 'se', 'pval', 'chr', 'pos']
        for col in numeric_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')

        # 2. 删除空值行
        existing_check_cols = [c for c in CRITICAL_COLS if c in df.columns]
        df.dropna(subset=existing_check_cols, inplace=True)
        
        # 3. 【关键修复】强制 chr 和 pos 转为整数 (Integer)
        # 这步是为了防止 1 变成 1.0，导致 "chr1.0" 这种错误
        if 'pos' in df.columns:
            df['pos'] = df['pos'].astype(int)
            # 备份原始坐标
            df['pos_original'] = df['pos'] 
            
        if 'chr' in df.columns:
            # 先转 int 去掉小数点，再转 string
            df['chr'] = df['chr'].astype(int) 
            df['temp_chr'] = df['chr'].astype(str) # 此时肯定是 "1", "2" 而不是 "1.0"
        
        cleaned_len = len(df)
        if cleaned_len < raw_len:
            print(f"   [清洗结果] 已删除 {raw_len - cleaned_len} 行坏数据。")
        
        if cleaned_len == 0:
            print("   [错误]清洗后数据为空，跳过此文件。")
            continue

        # ==================================================
        # 步骤 1: 统一转换为 hg19 (用于过滤)
        # ==================================================
        df_hg19 = None 
        
        if INPUT_BUILD == 'hg38':
            print(">> 步骤1 (判别:hg38): 执行转换 hg38 -> hg19 ...")
            new_pos_list = []
            valid_indices = []

            temp_chr_array = df['temp_chr'].values
            pos_array = df['pos'].values 
            indexes = df.index.values
            
            for idx, chrom, pos in tqdm(zip(indexes, temp_chr_array, pos_array), total=len(df), desc="LiftOver 38->19"):
                convert_chrom = f"chr{chrom}"
                new_coords = lo_38_to_19.convert_coordinate(convert_chrom, pos)
                if new_coords:
                    new_pos_list.append(new_coords[0][1])
                    valid_indices.append(idx)
            
            df_hg19 = df.loc[valid_indices].copy()
            df_hg19['pos'] = new_pos_list 
            print(f"   转换成功保留行数: {len(df_hg19)}")
            
        elif INPUT_BUILD == 'hg19':
            print(">> 步骤1 (判别:hg19): 输入已是 hg19，保留原坐标。")
            df_hg19 = df.copy()
        
        else:
            print(f"错误: INPUT_BUILD 配置错误")
            continue

        # ==================================================
        # 步骤 2 & 3: MHC 和 MAF 过滤
        # ==================================================
        print(">> 步骤2 & 3: 执行 MHC 和 MAF 过滤...")
        
        # 此时 df_hg19['temp_chr'] 已经是干净的 "1", "6" 等
        mask_mhc = (df_hg19['temp_chr'] == MHC_CHROM) & \
                   (df_hg19['pos'] >= MHC_START) & \
                   (df_hg19['pos'] <= MHC_END)
        df_hg19 = df_hg19[~mask_mhc]
        
        # MAF 过滤
        eaf_values = df_hg19['eaf'].values
        maf_values = np.minimum(eaf_values, 1 - eaf_values)
        mask_low_maf = maf_values < MAF_THRESHOLD
        df_hg19 = df_hg19[~mask_low_maf]
        
        print(f"   过滤后剩余行数: {len(df_hg19)}")

        # 删除临时列
        if 'temp_chr' in df_hg19.columns:
            del df_hg19['temp_chr']

        # ==================================================
        # 步骤 4: 输出到 1.hg19
        # ==================================================
        df_out_hg19 = df_hg19.drop(columns=['pos_original'], errors='ignore')
        out_path_hg19 = os.path.join("1.hg19", file_name)
        print(f">> 步骤4: 保存 hg19 文件 -> {out_path_hg19}")
        df_out_hg19.to_csv(out_path_hg19, sep='\t', index=False)

        # ==================================================
        # 步骤 5: 生成 FUMA 文件
        # ==================================================
        out_path_fuma = os.path.join("3.FUMA", file_name + ".gz")
        print(f">> 步骤5: 生成 FUMA 文件 -> {out_path_fuma}")
        
        fuma_cols = {
            'SNP': 'SNP', 'chr': 'CHR', 'pos': 'BP',
            'effect_allele': 'A1', 'other_allele': 'A2',
            'pval': 'P', 'beta': 'BETA', 'se': 'SE', 'samplesize': 'N'
        }
        missing_cols = [c for c in fuma_cols.keys() if c not in df_hg19.columns]
        if not missing_cols:
            df_fuma = df_hg19[list(fuma_cols.keys())].copy()
            df_fuma.rename(columns=fuma_cols, inplace=True)
            df_fuma.to_csv(out_path_fuma, sep='\t', index=False, compression='gzip')

        # ==================================================
        # 步骤 6: 生成 2.hg38 (修复版)
        # ==================================================
        print(">> 步骤6: 生成 hg38 版本数据 ...")
        out_path_hg38 = os.path.join("2.hg38", file_name)
        
        df_hg38_final = df_hg19.copy()
        
        if INPUT_BUILD == 'hg38':
            print("   (策略: 恢复原始 hg38 坐标)")
            df_hg38_final['pos'] = df_hg38_final['pos_original']
            if 'pos_original' in df_hg38_final.columns:
                del df_hg38_final['pos_original']
            df_hg38_final.to_csv(out_path_hg38, sep='\t', index=False)
            print(f"   保存成功 -> {out_path_hg38}")
            
        elif INPUT_BUILD == 'hg19':
            print("   (策略: 计算 hg19 -> hg38)")
            if 'pos_original' in df_hg38_final.columns:
                del df_hg38_final['pos_original']
            
            # 【重要】确保这里也是干净的整数转字符串
            df_hg38_final['chr'] = df_hg38_final['chr'].astype(int)
            df_hg38_final['temp_chr'] = df_hg38_final['chr'].astype(str)
            
            new_pos_38_list = []
            valid_indices_38 = []
            
            temp_chr_array = df_hg38_final['temp_chr'].values
            pos_array = df_hg38_final['pos'].values 
            indexes = df_hg38_final.index.values

            for idx, chrom, pos in tqdm(zip(indexes, temp_chr_array, pos_array), total=len(df_hg38_final), desc="LiftOver 19->38"):
                convert_chrom = f"chr{chrom}" # 这里现在肯定是 chr1, chr2 (不会是 chr1.0)
                new_coords = lo_19_to_38.convert_coordinate(convert_chrom, pos)
                if new_coords:
                    new_pos_38_list.append(new_coords[0][1])
                    valid_indices_38.append(idx)
            
            df_hg38_final = df_hg38_final.loc[valid_indices_38].copy()
            df_hg38_final['pos'] = new_pos_38_list
            
            if 'temp_chr' in df_hg38_final.columns:
                del df_hg38_final['temp_chr']
                
            df_hg38_final.to_csv(out_path_hg38, sep='\t', index=False)
            print(f"   保存成功 -> {out_path_hg38}")

    print("\n所有文件处理完成！")

if __name__ == "__main__":
    main()
