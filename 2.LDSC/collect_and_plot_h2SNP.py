import os
import glob
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# ================= 配置区域 =================
# 1. 输入：存放 log 文件的文件夹
INPUT_LOG_DIR = "2.h2"

# 2. 输出 1：基础数据汇总表 (CSV)
OUTPUT_CSV = "h2SNP.csv"

# 3. 输出 2：分组后的详细列表 (CSV)
OUTPUT_LIST = "h2SNP_grouped_list.csv"

# 4. 输出 3：透视矩阵表 (CSV, 适合论文直接使用)
OUTPUT_MATRIX = "h2SNP_matrix.csv"

# 5. 输出 4：分组森林图 (PDF)
OUTPUT_PLOT = "h2SNP_grouped_plot.pdf"
# ===========================================

def extract_from_log(file_path):
    """从单个 log 文件中提取关键指标"""
    data = {
        "File": os.path.basename(file_path).replace("_h2.log", ""),
        "h2": None, "h2_se": None, "Z_score": 0,
        "Lambda_GC": None, "Intercept": None, "Intercept_se": None
    }
    try:
        with open(file_path, 'r') as f:
            content = f.read()
            
            # 提取 h2 和 SE
            h2_match = re.search(r"Total Observed scale h2:\s+([-\d\.]+)\s+\(([\d\.]+)\)", content)
            if h2_match:
                data["h2"] = float(h2_match.group(1))
                data["h2_se"] = float(h2_match.group(2))
                if data["h2_se"] > 0:
                    data["Z_score"] = round(data["h2"] / data["h2_se"], 4)

            # 提取 Lambda GC
            lambda_match = re.search(r"Lambda GC:\s+([-\d\.]+)", content)
            if lambda_match:
                data["Lambda_GC"] = float(lambda_match.group(1))

            # 提取 Intercept
            int_match = re.search(r"Intercept:\s+([-\d\.]+)\s+\(([\d\.]+)\)", content)
            if int_match:
                data["Intercept"] = float(int_match.group(1))
                data["Intercept_se"] = float(int_match.group(2))
                
    except Exception as e:
        print(f"读取错误 {file_path}: {e}")
    return data

def parse_filename(name):
    """
    解析文件名逻辑：Source_Phenotype
    例如: UKB_Height -> Source: UKB, Phenotype: Height
    使用第一个下划线进行分割 (从左往右切 1 刀)
    """
    if '_' in name:
        parts = name.split('_', 1) # 这里根据需要调整：split(左切) 或 rsplit(右切)
        return parts[0], parts[1]
    return "Unknown", name

def main():
    # ==================================================
    # 第一步：数据提取 (Extract)
    # ==================================================
    print(f"1. 正在从 {INPUT_LOG_DIR} 提取数据...")
    if not os.path.exists(INPUT_LOG_DIR):
        print(f"错误: 找不到目录 {INPUT_LOG_DIR}")
        return

    log_files = glob.glob(os.path.join(INPUT_LOG_DIR, "*.log"))
    all_results = []
    
    for log_file in log_files:
        res = extract_from_log(log_file)
        if res["h2"] is not None:
            all_results.append(res)

    if not all_results:
        print("未提取到有效数据，程序终止。")
        return

    # 生成基础 DataFrame
    df = pd.DataFrame(all_results)
    
    # 增加分组列
    df[['Source', 'Phenotype']] = df['File'].apply(lambda x: pd.Series(parse_filename(x)))
    
    # 排序：按表型字母顺序，再按来源顺序
    df = df.sort_values(by=['Phenotype', 'Source'])
    
    # 保存基础 CSV
    cols_order = ['File', 'Source', 'Phenotype', 'h2', 'h2_se', 'Z_score', 'Intercept', 'Intercept_se', 'Lambda_GC']
    df[cols_order].to_csv(OUTPUT_CSV, index=False)
    print(f"✔ 原始数据提取完成: {OUTPUT_CSV}")

    # ==================================================
    # 第二步：生成矩阵表 (Pivot Table)
    # ==================================================
    print("2. 正在生成矩阵表...")
    # 创建显示列 "h2 (se)"
    df['Display'] = df.apply(lambda row: f"{row['h2']:.4f} ({row['h2_se']:.4f})", axis=1)
    
    try:
        pivot_df = df.pivot(index='Phenotype', columns='Source', values='Display')
        pivot_df.to_csv(OUTPUT_MATRIX)
        print(f"✔ 矩阵表已生成: {OUTPUT_MATRIX} (适合论文对比)")
        
        # 保存一份带分组的详细列表
        df[cols_order].to_csv(OUTPUT_LIST, index=False)
        print(f"✔ 分组列表已生成: {OUTPUT_LIST}")
    except Exception as e:
        print(f"⚠️ 矩阵表生成失败: {e}")

    # ==================================================
    # 第三步：分组绘图 (Visualization)
    # ==================================================
    print("3. 正在生成分组森林图...")
    
    phenotypes = df['Phenotype'].unique()
    sources = df['Source'].unique()
    
    # 动态计算颜色
    colors = sns.color_palette("husl", len(sources))
    color_map = dict(zip(sources, colors))

    # 动态计算图片高度
    fig_height = max(6, len(phenotypes) * len(sources) * 0.5)
    plt.figure(figsize=(10, fig_height))

    # 绘图参数
    y_base = np.arange(len(phenotypes)) 
    bar_height = 0.6  
    step = bar_height / len(sources)

    # 循环绘制
    for i, pheno in enumerate(phenotypes):
        group_data = df[df['Phenotype'] == pheno]
        
        for j, src in enumerate(sources):
            row = group_data[group_data['Source'] == src]
            if not row.empty:
                # 计算 Y 轴偏移量
                offset = (j - len(sources)/2) * step + (step/2)
                y_pos = i + offset
                
                h2_val = row['h2'].values[0]
                se_val = row['h2_se'].values[0]
                
                # 画点和误差棒
                plt.errorbar(h2_val, y_pos, xerr=1.96*se_val, 
                             fmt='o', color=color_map[src], 
                             ecolor='gray', capsize=3, markersize=6,
                             label=src if i == 0 else "") # 防止图例重复

    # 设置轴和标签
    plt.yticks(y_base, phenotypes, fontsize=11, fontweight='bold')
    plt.xlabel('SNP Heritability ($h^2$)', fontsize=12)
    plt.title('Heritability Estimates (Grouped by Phenotype)', fontsize=14)
    
    plt.axvline(x=0, color='black', linestyle='--', linewidth=0.8)
    plt.legend(title="Source (Prefix)", loc='upper right')
    plt.grid(axis='x', linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, format='pdf')
    print(f"✔ 绘图完成: {OUTPUT_PLOT}")
    print("="*40)
    print("所有任务全部完成！")

if __name__ == "__main__":
    main()
