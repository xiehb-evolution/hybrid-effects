#!/bin/bash

source activate vep_env
export LD_LIBRARY_PATH=/data_group/xiehaibing/xiehaibing1/f2/htslib-1.20:$LD_LIBRARY_PATH

process_bam() {
  local id=$1
  local bam_dir=$2
  local output_dir=$3
  local bam_file="${bam_dir}/${id}.bam"
  local name_bam="${output_dir}/${id}.name.bam"
  local pass_bam="${output_dir}/${id}.pass.bam"
  local pass_sorted_bam="${output_dir}/${id}.pass.sorted.bam"
  local fail_sorted_bam="${output_dir}/${id}.fail.sorted.bam"
  local log_file="${output_dir}/${id}.log"

  # 检查 BAM 文件是否存在
  if [[ ! -f "$bam_file" ]]; then
    echo "Error: ${bam_file} does not exist!"
    return 1
  fi

  # 检查目标文件是否已存在，如果存在则跳过处理
  if [[ -f "$pass_sorted_bam" && -f "$fail_sorted_bam" ]]; then
    echo "Output files for ${id} already exist, skipping."
    return 0
  fi

  # 检查中间文件是否存在，如果存在则认为其他进程正在处理
  if [[ -f "$name_bam" || -f "$pass_bam" ]]; then
    echo "Another process is already handling ${id}, skipping."
    return 0
  fi

  # 检查输出目录是否存在，不存在则创建
  if [[ ! -d "$output_dir" ]]; then
    echo "Output directory ${output_dir} does not exist. Creating it."
    mkdir -p "$output_dir"
  fi

  # 开始处理 BAM 文件
  echo "Processing ${id}..."
  samtools sort -n "${bam_file}" -o "${name_bam}"
  ./filter_bam "${name_bam}" "${pass_bam}" "${fail_bam}" 20 1000 15 2 | gzip -c - > "${log_file}"
  samtools sort "${pass_bam}" -o "${pass_sorted_bam}" --write-index
  samtools sort "${fail_bam}" -o "${fail_sorted_bam}" --write-index

  # 清理中间文件
  rm "${name_bam}" "${pass_bam}" "${fail_bam}"
}

# 检查是否提供了必要的参数
if [[ -z "$1" || -z "$2" || -z "$3" ]]; then
  echo "Error: Missing arguments."
  echo "Usage: $0 <id_list_file> <bam_directory> <output_directory>"
  exit 1
fi

id_list_file=$1
bam_dir=$2
output_dir=$3

# 检查 id_list_file 是否存在
if [[ ! -f "$id_list_file" ]]; then
  echo "Error: ID list file ${id_list_file} does not exist!"
  exit 1
fi

# 逐行读取 id_list_file，并调用 process_bam 函数并行处理
while IFS= read -r id; do
  if [[ -n "$id" ]]; then  # 检查是否为空行
    process_bam "$id" "$bam_dir" "$output_dir"
  fi
done < "$id_list_file"
