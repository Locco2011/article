#!/bin/bash

# ================= 配置区域 =================
# 1. 输入文件的路径
INPUT_DIR="../1.QC/1.hg19"

# 2. 输出文件的文件夹名称
OUTPUT_DIR="1.sumstats"

# 3. munge_sumstats.py 的绝对路径
MUNGE_SCRIPT="/home/chenguanglei2011/ldsc/munge_sumstats.py"

# 4. w_hm3.snplist 的绝对路径
HM3_PATH="/home/chenguanglei2011/ldsc/eur_w_ld_chr/w_hm3.snplist"

# 5. N 列的列名
N_COL="samplesize"
# ===========================================

# 检查输入目录是否存在
if [ ! -d "$INPUT_DIR" ]; then
    echo "错误: 找不到输入目录 $INPUT_DIR"
    exit 1
fi

# 检查脚本是否存在
if [ ! -f "$MUNGE_SCRIPT" ]; then
    echo "错误: 找不到脚本文件 $MUNGE_SCRIPT"
    exit 1
fi

# 如果输出目录不存在，则创建它
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "创建输出目录: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
fi

echo "开始批量处理..."
echo "输出位置: $OUTPUT_DIR/"

# 循环读取目录下的所有 .txt 文件
for file in "$INPUT_DIR"/*.txt; do
    # 检查是否真的有文件
    [ -e "$file" ] || continue

    # 1. 提取文件名 (例如: Fin_ABC.txt)
    filename=$(basename "$file")
    
    # 2. 提取前缀 (例如: Fin_ABC)
    prefix="${filename%.*}"
    
    # 3. 定义预期的输出完整路径
    # 注意: munge_sumstats 会自动添加 .sumstats.gz，所以我们只需要指定目录/前缀
    expected_output="$OUTPUT_DIR/${prefix}.sumstats.gz"

    # ==================================================
    # 4. 断点续传 (优化版)
    # ==================================================
    # 使用 -s 参数：检查文件是否存在 且 文件大小大于0
    # 这样可以避免上次运行中断产生空文件，导致本次被错误跳过的情况
    if [ -s "$expected_output" ]; then
        echo ">>> [跳过] 文件已存在且完整: $expected_output"
        continue
    fi
    # ==================================================

    echo "=================================================="
    echo "正在处理: $prefix"
    
    # 5. 执行命令
    # 注意 --out 参数现在包含了输出目录
    python "$MUNGE_SCRIPT" \
        --sumstats "$file" \
        --N-col "$N_COL" \
        --out "$OUTPUT_DIR/$prefix" \
        --merge-alleles "$HM3_PATH"

    # 检查结果
    # 再次检查生成的文件是否存在且不为空
    if [ -s "$expected_output" ]; then
        echo "✔ 成功: $expected_output"
    else
        echo "✘ 失败: $filename (可能未生成文件或文件为空)"
    fi

done

echo "所有任务完成！"
