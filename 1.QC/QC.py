import pandas as pd
from pyliftover import LiftOver
from tqdm import tqdm
import os
import glob
import numpy as np

# ================= 配置区域 =================
# 1. 输入数据的文件夹路径
INPUT_FOLDER = "./0.Raw/Fin" 

# 2. 定义 MHC 区域 (基于 hg19 坐标)
MHC_CHROM = '6'  
MHC_START = 28477797
MHC_END = 33448354

# 3. MAF 阈值
MAF_THRESHOLD = 0.01

# 4. 需要检查空值的关键列
CRITICAL_COLS = ['SNP', 'effect_allele', 'other_allele', 'eaf', 'beta', 'se', 'pval', 'chr', 'pos']
# ===========================================

def main():
    # --- 0. 创建输出目录 ---
    dirs = ["1.hg19", "2.hg38", "3.FUMA"]
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    # --- 1. 初始化 LiftOver ---
    print("正在初始化转换链...")
    try:
        lo_38_to_19 = LiftOver('hg38', 'hg19')
        lo_19_to_38 = LiftOver('hg19', 'hg38')
    except Exception as e:
        print(f"初始化 LiftOver 失败: {e}")
        return

    # 获取所有 txt 文件
    files = glob.glob(os.path.join(INPUT_FOLDER, "*.txt"))
    print(f"共发现 {len(files)} 个文件需要处理。")

    for file_path in files:
        file_name = os.path.basename(file_path)
        print(f"\n{'='*50}")
        print(f"正在处理文件: {file_name}")
        
        # --- 读取数据 ---
        try:
            # 这里的 sep=r'\s+' 能很好地处理空格或制表符
            df = pd.read_csv(file_path, sep=r'\s+', low_memory=False)
        except Exception as e:
            print(f"读取失败，跳过: {e}")
            continue

        raw_len = len(df)
        
        # ==================================================
        # 【新增步骤】 深度数据清洗
        # ==================================================
        print(f">> 数据清洗前行数: {raw_len}")
        
        # 1. 确保数值列真的是数值 (处理列错位导致的字符串问题)
        # 如果某一行发生粘连(如 0.01230.456)，强制转为 NaN
        numeric_cols = ['eaf', 'beta', 'se', 'pval', 'chr', 'pos']
        for col in numeric_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')

        # 2. 检查关键列是否存在，如果存在则进行空值删除
        existing_check_cols = [c for c in CRITICAL_COLS if c in df.columns]
        
        if len(existing_check_cols) < len(CRITICAL_COLS):
             print(f"   [警告] 文件缺少部分关键列，仅检查: {existing_check_cols}")
        
        # 3. 删除包含 NaN 的行
        df.dropna(subset=existing_check_cols, inplace=True)
        
        # 4. 确保 pos 和 chr 转为整数 (前面转numeric后可能是float)
        if 'pos' in df.columns:
            df['pos'] = df['pos'].astype(int)
        if 'chr' in df.columns:
            # chr有时可能是 X, Y，转换 int 可能会报错，这里假设是常染色体或者转字符串处理
            # 为了保险，我们先统一转字符串去掉 'chr'，再尝试转 int (如果全是数字)
            df['temp_chr'] = df['chr'].astype(str).str.replace('chr', '', regex=False)
        
        cleaned_len = len(df)
        if cleaned_len < raw_len:
            print(f"   [清洗结果] 已删除 {raw_len - cleaned_len} 行坏数据 (包含空值或格式错误)。")
        
        if cleaned_len == 0:
            print("   [错误]清洗后数据为空，跳过此文件。")
            continue

        # ==================================================
        # 步骤 1: hg38 -> hg19 转换 (使用 zip 提速)
        # ==================================================
        print(">> 步骤1: 转换坐标 hg38 -> hg19 ...")
        new_pos_list = []
        valid_indices = []

        temp_chr_array = df['temp_chr'].values
        pos_array = df['pos'].values
        indexes = df.index.values
        
        # 使用 zip 遍历
        for idx, chrom, pos in tqdm(zip(indexes, temp_chr_array, pos_array), total=len(df), desc="LiftOver"):
            convert_chrom = f"chr{chrom}"
            new_coords = lo_38_to_19.convert_coordinate(convert_chrom, pos)
            
            if new_coords:
                new_pos_list.append(new_coords[0][1])
                valid_indices.append(idx)
        
        df_hg19 = df.loc[valid_indices].copy()
        df_hg19['pos'] = new_pos_list 
        
        # ==================================================
        # 步骤 2: 删除 MHC 区域 (chr6: 28,477,797–33,448,354)
        # ==================================================
        print(">> 步骤2: 过滤 MHC 区域...")
        df_hg19['temp_chr'] = df_hg19['chr'].astype(str).str.replace('chr', '', regex=False)
        
        mask_mhc = (df_hg19['temp_chr'] == MHC_CHROM) & \
                   (df_hg19['pos'] >= MHC_START) & \
                   (df_hg19['pos'] <= MHC_END)
        
        n_mhc = mask_mhc.sum()
        df_hg19 = df_hg19[~mask_mhc]
        print(f"   已移除 MHC 区域变异数: {n_mhc}")

        # ==================================================
        # 步骤 3: 计算 MAF 并删除 < 0.01
        # ==================================================
        print(">> 步骤3: MAF 过滤 (< 0.01)...")
        # eaf 已经是 numeric 且非空了
        eaf_values = df_hg19['eaf'].values
        maf_values = np.minimum(eaf_values, 1 - eaf_values)
        
        mask_low_maf = maf_values < MAF_THRESHOLD
        n_low_maf = mask_low_maf.sum()
        
        df_hg19 = df_hg19[~mask_low_maf]
        print(f"   已移除 MAF < 0.01 变异数: {n_low_maf}")
        print(f"   清洗后剩余总行数: {len(df_hg19)}")

        # 删除临时列
        if 'temp_chr' in df_hg19.columns:
            del df_hg19['temp_chr']

        # ==================================================
        # 步骤 4: 输出到 1.hg19 文件夹
        # ==================================================
        out_path_hg19 = os.path.join("1.hg19", file_name)
        print(f">> 步骤4: 保存 hg19 文件 -> {out_path_hg19}")
        df_hg19.to_csv(out_path_hg19, sep='\t', index=False)

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
        else:
            print(f"   [警告] 缺少列 {missing_cols}，跳过 FUMA 文件生成。")

        # ==================================================
        # 步骤 6: 将最终结果转回 hg38 -> 2.hg38
        # ==================================================
        print(">> 步骤6: 将最终结果转回 hg38 (hg19 -> hg38) ...")
        
        df_hg38_final = df_hg19.copy()
        df_hg38_final['temp_chr'] = df_hg38_final['chr'].astype(str).str.replace('chr', '', regex=False)
        
        new_pos_38_list = []
        valid_indices_38 = []
        
        temp_chr_array_38 = df_hg38_final['temp_chr'].values
        pos_array_38 = df_hg38_final['pos'].values 
        indexes_38 = df_hg38_final.index.values

        for idx, chrom, pos in tqdm(zip(indexes_38, temp_chr_array_38, pos_array_38), total=len(df_hg38_final), desc="LiftOver 19->38"):
            convert_chrom = f"chr{chrom}"
            new_coords = lo_19_to_38.convert_coordinate(convert_chrom, pos)
            
            if new_coords:
                new_pos_38_list.append(new_coords[0][1])
                valid_indices_38.append(idx)
        
        df_hg38_final = df_hg38_final.loc[valid_indices_38].copy()
        df_hg38_final['pos'] = new_pos_38_list
        
        if 'temp_chr' in df_hg38_final.columns:
            del df_hg38_final['temp_chr']
            
        out_path_hg38 = os.path.join("2.hg38", file_name)
        print(f"   保存 hg38 文件 -> {out_path_hg38}")
        df_hg38_final.to_csv(out_path_hg38, sep='\t', index=False)

    print("\n所有文件处理完成！")

if __name__ == "__main__":
    main()
