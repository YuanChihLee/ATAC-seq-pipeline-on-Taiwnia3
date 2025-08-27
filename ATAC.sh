#!/bin/bash
#SBATCH -A you_account      # Specify the project account
#SBATCH -J ATACseq_mm39_Final # Specify the job name
#SBATCH -p ct56          # Specify the partition
#SBATCH -n 14             # Specify the total number of cores
#SBATCH --cpus-per-task=4 # Specify the number of CPUs per task
#SBATCH -o %j.out         # Standard output file
#SBATCH -e %j.err         # Standard error file

# Stop the script immediately if any command fails or an unset variable is used
set -e
set -u
set -o pipefail

# === Global Parameters ===
# --- User-defined paths ---
proj_root="/work/you_account/ATAC-seq" # Your ATAC-seq project root directory
pkg_dir="/work/you_account/CUTandTAG/CUTandTag_pkg" # Directory for reference data and packages
raw_dir="$proj_root/ATAC"       # Directory for raw fastq files
filelist="$proj_root/filelist.txt"   # Sample list file

# --- Reference Genome and Blacklist Files (mm39) ---
GENOME_INDEX="$pkg_dir/index/mouse/mm39"
GENOME_SIZE="mm"
CHROM_SIZES="$pkg_dir/mm39.chrom.sizes"
BLACKLIST_FILE="$pkg_dir/mm39.excluderanges_2.bed"

# --- Automatically generated directories ---
trim_dir="$proj_root/01_trimmed"
map_dir="$proj_root/02_mapping"
final_bam_dir="$proj_root/03_final_bam"
peak_dir="$proj_root/04_peaks"
bw_dir="$proj_root/05_bigwig"
log_dir="$proj_root/logs"
qc_dir="$proj_root/QC_reports"
tmp_dir="$proj_root/00_tmp"

# Explicitly specify the path to your local MACS3 installation
export MACS3_PATH="$HOME/.local/bin/macs3"

# Export global variables for use with GNU Parallel
export proj_root pkg_dir filelist GENOME_INDEX GENOME_SIZE CHROM_SIZES BLACKLIST_FILE \
       raw_dir trim_dir map_dir final_bam_dir peak_dir bw_dir log_dir qc_dir tmp_dir MACS3_PATH

# Create required directories
mkdir -p "$trim_dir" "$map_dir" "$tmp_dir" "$final_bam_dir" "$peak_dir" "$bw_dir" "$log_dir" "$qc_dir"

# Clear the previous QC summary file
> "$qc_dir/QC_summary.tsv"

# === Load Environment Modules ===
init_env() {
  source /etc/profile.d/modules.sh
  module purge
  module load old-module
  module load biology/Trimmomatic/0.39
  module load biology/SAMTOOLS/1.18
  module load biology/bowtie2/2.4.2
  module load biology/BEDTOOLS/2.31.1
  module load biology/Picard/2.27.4
  export PATH=$PATH:/opt/ohpc/Taiwania3/pkg/biology/UCSC_Utilities/Utilities_20180515
}
export -f init_env

# === [1. Trim Reads - Trimming] ===
run_trim() {
  init_env
  local sample_id="$1"
  local i5="$2"
  local i7="$3"
  local fq1_candidates=("$raw_dir/${sample_id}"*R1*.fastq.gz)
  local fq2_candidates=("$raw_dir/${sample_id}"*R2*.fastq.gz)
  if [[ ${#fq1_candidates[@]} -ne 1 ]]; then echo "[ERROR] Found ${#fq1_candidates[@]} R1 files for ${sample_id}." >&2; return 1; fi
  if [[ ${#fq2_candidates[@]} -ne 1 ]]; then echo "[ERROR] Found ${#fq2_candidates[@]} R2 files for ${sample_id}." >&2; return 1; fi
  local fq1="${fq1_candidates[0]}"
  local fq2="${fq2_candidates[0]}"
  local adaptor_file="$pkg_dir/v2_Ad1_${i5}/v2_Ad2_${i7}.fa"
  if [[ ! -f "$adaptor_file" ]]; then echo "[ERROR] Adaptor file $adaptor_file not found for sample ${sample_id}." >&2; return 1; fi
  
  echo "[INFO] [${sample_id}] Starting Trimmomatic..."
  trimmomatic PE -threads 4 -phred33 "$fq1" "$fq2" \
    "$trim_dir/${sample_id}_R1.paired.fastq.gz" "$trim_dir/${sample_id}_R1.unpaired.fastq.gz" \
    "$trim_dir/${sample_id}_R2.paired.fastq.gz" "$trim_dir/${sample_id}_R2.unpaired.fastq.gz" \
    ILLUMINACLIP:"$adaptor_file":2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:30 \
    2> "$log_dir/${sample_id}.trimmomatic.log"
  echo "[INFO] [${sample_id}] Trimmomatic completed."
}
export -f run_trim

# === [2. Sequence Alignment] ===
run_alignment() {
  init_env
  local sample_id="$1"
  local fq1="$trim_dir/${sample_id}_R1.paired.fastq.gz"
  local fq2="$trim_dir/${sample_id}_R2.paired.fastq.gz"
  local sam_file="$map_dir/${sample_id}.sam"
  local log_file="$log_dir/${sample_id}.bowtie2.log"

  if [[ ! -f "$fq1" || ! -f "$fq2" ]]; then echo "[WARN] Missing trimmed FASTQ for $sample_id. Skipping alignment." >&2; return 1; fi

  echo "[INFO] [${sample_id}] Aligning to mm39 genome..."
  bowtie2 --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 2000 -p 4 \
    -x "$GENOME_INDEX" -1 "$fq1" -2 "$fq2" \
    -S "$sam_file" 2> "$log_file"
  
  echo "[INFO] [${sample_id}] Bowtie2 alignment finished. Checking for SAM file..."
  if [ -s "$sam_file" ]; then
      echo "[SUCCESS] [${sample_id}] SAM file '$sam_file' created successfully and is not empty."
  else
      echo "[ERROR] [${sample_id}] SAM file was NOT created or is empty. Alignment may have failed. Please check log: '$log_file'." >&2
      return 1
  fi
  
  local aln_rate=$(grep "overall alignment rate" "$log_file" | cut -d' ' -f1)
  echo -e "${sample_id}\tAlignmentRate\t${aln_rate}" > "$log_dir/${sample_id}.qc_temp.txt"
}
export -f run_alignment

# === [3. Post-alignment Processing and Filtering] ===
run_postprocess() {
    init_env
    local sample_id="$1"
    local sam_file="$map_dir/${sample_id}.sam"
    local qc_temp_file="$log_dir/${sample_id}.qc_temp.txt"
    
    local sorted_bam="$map_dir/${sample_id}.sorted.bam"
    local rmdup_bam="$map_dir/${sample_id}.rmDup.bam"
    local rmdup_metrics="$qc_dir/${sample_id}.rmDup.metrics.txt"
    local temp_filtered_bam="${tmp_dir}/${sample_id}.temp_filtered.bam"
    local final_bam="$final_bam_dir/${sample_id}.final.bam"

    echo "[INFO] [${sample_id}] Starting post-alignment processing..."
    samtools view -bS "$sam_file" | samtools sort -@ 4 -T "$tmp_dir/${sample_id}.sort" -o "$sorted_bam"
    
    # As requested, do not remove the SAM file.
    # rm "$sam_file"
    
    echo "[INFO] [${sample_id}] Removing PCR duplicates with Picard..."
    picard MarkDuplicates I="$sorted_bam" O="$rmdup_bam" \
        REMOVE_DUPLICATES=true METRICS_FILE="$rmdup_metrics" ASSUME_SORT_ORDER=coordinate

    echo "[INFO] [${sample_id}] Indexing de-duplicated BAM..."
    samtools index "$rmdup_bam"

    if [ -f "$rmdup_metrics" ]; then
        local dup_rate=$(grep -A1 "LIBRARY" "$rmdup_metrics" | tail -n1 | cut -f9)
        echo -e "${sample_id}\tDuplicationRate\t${dup_rate}" >> "$qc_temp_file"
    fi

    echo "[INFO] [${sample_id}] Filtering reads..."
    local total_reads=$(samtools view -c "$rmdup_bam")
    local chrM_reads=$(samtools view -c "$rmdup_bam" chrM MT 2>/dev/null || echo 0)
    if [ "$total_reads" -gt 0 ]; then
        local chrM_rate=$(awk "BEGIN {print ($chrM_reads / $total_reads) * 100}")
        echo -e "${sample_id}\tMitoContaminationRate\t${chrM_rate}%" >> "$qc_temp_file"
    fi
    
    echo "[INFO] [${sample_id}] Applying filters (MAPQ >= 10, proper pair, no unmapped/secondary)..."
    samtools view -b -F 1804 -f 2 -q 10 "$rmdup_bam" > "$temp_filtered_bam"
    
    echo "[INFO] [${sample_id}] Indexing temporary filtered BAM..."
    samtools index "$temp_filtered_bam" # This is a critical fix

    echo "[INFO] [${sample_id}] Removing mitochondrial and blacklisted reads..."
    local non_mito_chroms=$(samtools idxstats "$temp_filtered_bam" | cut -f1 | grep -v -E "^(chrM|MT|\*)$")
    samtools view -b "$temp_filtered_bam" ${non_mito_chroms} | \
        bedtools intersect -v -a - -b "$BLACKLIST_FILE" > "$final_bam"
    
    if [ ! -s "$final_bam" ]; then
        echo "[ERROR] [${sample_id}] Final BAM file is empty after all filters." >&2
        local temp_count=$(samtools view -c "$temp_filtered_bam")
        echo "[DEBUG] [${sample_id}] Reads before blacklist/mito filtering: ${temp_count}" >&2
        return 1
    fi
    
    samtools index "$final_bam"

    echo "[INFO] [${sample_id}] Cleaning up intermediate files..."
    rm "$sorted_bam" "$rmdup_bam" "${rmdup_bam}.bai" "$temp_filtered_bam" "${temp_filtered_bam}.bai"
    echo "[INFO] [${sample_id}] Post-processing completed. Final BAM is at ${final_bam}"
}
export -f run_postprocess

# === [4. Peak Calling & BigWig Generation] ===
run_peak_and_bw() {
    init_env
    local sample_id="$1"
    local final_bam="$final_bam_dir/${sample_id}.final.bam"
    local qc_temp_file="$log_dir/${sample_id}.qc_temp.txt"

    if [ ! -s "$final_bam" ]; then
        echo "[ERROR] [${sample_id}] Final BAM is missing or empty. Skipping peak calling and BigWig generation." >&2
        return 1
    fi

    echo "[INFO] [${sample_id}] Generating BigWig file..."
    local bedgraph="$bw_dir/${sample_id}.bedgraph"
    local final_bw="$bw_dir/${sample_id}.bw"
    
    local total_reads_in_final_bam=$(samtools view -c "$final_bam")
    if [ ! "$total_reads_in_final_bam" -gt 0 ]; then
        echo "[WARN] [${sample_id}] No reads in final BAM. Cannot generate BigWig." >&2
    else
        local scale_factor=$(echo "1000000 / $total_reads_in_final_bam" | bc -l)
        bedtools genomecov -bg -ibam "$final_bam" -scale "$scale_factor" | sort -k1,1 -k2,2n > "$bedgraph"
        bedGraphToBigWig "$bedgraph" "$CHROM_SIZES" "$final_bw"
    fi
    
    echo "[INFO] [${sample_id}] Calling peaks with MACS3 for ATAC-seq..."
    "$MACS3_PATH" callpeak -t "$final_bam" -f BAMPE -n "${sample_id}" \
      -g "$GENOME_SIZE" --outdir "$peak_dir" \
      --nomodel --shift -100 --extsize 200 -q 0.01 --keep-dup all
    
    echo "[INFO] [${sample_id}] Calculating FrIP Score..."
    local peak_file="$peak_dir/${sample_id}_peaks.narrowPeak"
    if [ -f "$peak_file" ] && [ "$total_reads_in_final_bam" -gt 0 ]; then
        local reads_in_peaks=$(samtools view -c -L "$peak_file" "$final_bam")
        local frip_score=$(awk "BEGIN {print ($reads_in_peaks / $total_reads_in_final_bam) * 100}")
        echo -e "${sample_id}\tFrIP_Score\t${frip_score}%" >> "$qc_temp_file"
    else
        echo "[WARN] MACS3 peak file not found or no reads in BAM for ${sample_id}. Skipping FrIP calculation." >&2
    fi
    
    if [ -f "$bedgraph" ]; then
        rm "$bedgraph"
    fi
    echo "[INFO] [${sample_id}] Peak Calling & BigWig generation completed."
}
export -f run_peak_and_bw

# === Main Execution Logic ===
N_JOBS=${SLURM_NPROCS:-1}

echo "====== [Step 1/4] Starting Read Trimming ======"
#cat "$filelist" | parallel -j $N_JOBS --colsep '\s+' --joblog "$log_dir/parallel_trim.log" run_trim {1} {2} {3}

echo "====== [Step 2/4] Starting Sequence Alignment ======"
#cut -f1 "$filelist" | parallel -j $N_JOBS --joblog "$log_dir/parallel_align.log" run_alignment {}

echo "====== [Step 3/4] Starting Post-processing and Filtering ======"
cut -f1 "$filelist" | parallel -j $N_JOBS --joblog "$log_dir/parallel_postprocess.log" run_postprocess {}

echo "====== [Step 4/4] Starting Peak Calling and BigWig Generation ======"
cut -f1 "$filelist" | parallel -j $N_JOBS --joblog "$log_dir/parallel_peak.log" run_peak_and_bw {}

# === Final Step: Consolidate QC Reports ===
echo "====== Consolidating all QC reports ======"
{
  echo -e "SampleID\tMetric\tValue"
  find "$log_dir" -name "*.qc_temp.txt" -exec cat {} +
} > "$qc_dir/QC_summary.tsv"
find "$log_dir" -name "*.qc_temp.txt" -delete

echo "====== ATAC-seq pipeline finished successfully! ======"
