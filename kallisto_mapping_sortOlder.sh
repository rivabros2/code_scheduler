#!/usr/bin/bash
#SBATCH --job-name=kallisto
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=50
#SBATCH --time=1:00:00
#SBATCH --output=kallisto_mapping50cpu.%j.out
#SBATCH --error=kallisto_mapping.%j50cpu.err

# Directories
# FASTQ_DIR="/home/shared/epilab/BiomarkersML/fastqNCBI_SoyBean/fastq_temp_batch2"


### kallisto_mapping_sortOlder.sh  
###    squeue -u $USER --name=kallisto | tee -a $LOG_FILE


FASTQ_DIR="/home/shared/epilab/BiomarkersML/fastqNCBI_SoyBean/fastq_temp_batch2"
KALLISTO_OUTPUT_DIR="/home/shared/epilab/BiomarkersML/fastqNCBI_SoyBean/kallisto_results_batch4"
KALLISTO_INDEX="/home/shared/epilab/BiomarkersML/refGenomeSoyBean/GCF_000004515.6/87848rnas_Gmax_v4"

mkdir -p "$KALLISTO_OUTPUT_DIR"

# Loop through files sorted by modification time (oldest first)
for fastq1 in $(ls -1t --reverse "$FASTQ_DIR"/*_1.fastq); do
  fastq2="${fastq1/_1.fastq/_2.fastq}"
  base_name=$(basename "$fastq1" _1.fastq)
  
### Cuando el file contiene el prefijo SAM, significa que se bajo completamente y fue re-nombrado, asegura solo usar archivos completos
  if [[ "$fastq1" == *SAM* ]] && [[ -f "$fastq2" ]]; then
    echo "Mapping paired-end files for $base_name"
    /home/mrivarola/bin/kallisto quant -i "$KALLISTO_INDEX" --threads 50 \
      -o "$KALLISTO_OUTPUT_DIR/${base_name}" "$fastq1" "$fastq2" 2> "$KALLISTO_OUTPUT_DIR/${base_name}.err"
    if [ $? -eq 0 ]; then
      rm -f "$fastq1" "$fastq2"
    else
      echo "Mapping failed for $base_name. Check logs."
    fi
  elif [[ "$fastq1" == *SAM* ]]; then
    single_fastq="${fastq1/_1.fastq/.fastq}"
    if [ -f "$single_fastq" ]; then
      echo "Mapping single-end file for $base_name"
      /home/mrivarola/bin/kallisto quant -i "$KALLISTO_INDEX" --threads 50 --single -l 200 -s 20 \
        -o "$KALLISTO_OUTPUT_DIR/${base_name}" "$single_fastq" 2> "$KALLISTO_OUTPUT_DIR/${base_name}.err"
      if [ $? -eq 0 ]; then
        rm -f "$single_fastq"
      else
        echo "Mapping failed for $base_name. Check logs."
      fi
    fi
  else
    echo "Skipping file $fastq1 (not SAM or SAM)."
  fi
done


