#!/usr/bin/bash

# Input file containing run and sample names
INPUT_FILE="/home/shared/epilab/BiomarkersML/fastqNCBI_SoyBean/2240part_03.txt"

# Directory to store temporary FASTQ files
FASTQ_DIR="/home/shared/epilab/BiomarkersML/fastqNCBI_SoyBean/fastq_temp_batch2"

# Log directory
LOG_DIR="/home/shared/epilab/BiomarkersML/fastqNCBI_SoyBean/logs"
mkdir -p "$FASTQ_DIR" "$LOG_DIR"
exec > >(tee -a "$LOG_DIR/job_$SLURM_JOB_ID.log") 2>&1

# Retry settings
MAX_RETRIES=3

# Loop through each line in the input file
while IFS=',' read -r run_name sample_name; do
  echo "Processing run: $run_name for sample: $sample_name"

  # Attempt to download FASTQ files
  attempt=0
  success=0
  while [ $attempt -lt $MAX_RETRIES ]; do
    echo "Attempting fastq-dump for $run_name (Attempt: $((attempt + 1)))"
    /home/mrivarola/sratoolkit.3.1.1-ubuntu64/bin/fastq-dump --split-3 --outdir "$FASTQ_DIR" "$run_name" && success=1 && break
    attempt=$((attempt + 1))
    echo "fastq-dump failed for $run_name. Retrying..."
  done

  # Fallback to prefetch if fastq-dump fails
  if [ $success -eq 0 ]; then
    echo "Trying prefetch for $run_name..."
    /home/mrivarola/sratoolkit.3.1.1-ubuntu64/bin/prefetch "$run_name"
    if [ $? -eq 0 ]; then
      echo "prefetch succeeded for $run_name. Attempting fastq-dump again..."
      /home/mrivarola/sratoolkit.3.1.1-ubuntu64/bin/fastq-dump --split-3 --outdir "$FASTQ_DIR" "$run_name"
      if [ $? -ne 0 ]; then
        echo "Error downloading run $run_name even after prefetch." >&2
        continue
      fi
    else
      echo "prefetch failed for $run_name. Skipping..." >&2
      continue
    fi
  fi

  # Rename files for clarity
  fastq1="$FASTQ_DIR/${run_name}_1.fastq"
  fastq2="$FASTQ_DIR/${run_name}_2.fastq"
  single_fastq="$FASTQ_DIR/${run_name}.fastq"

  renamed_fastq1="$FASTQ_DIR/${sample_name}_${run_name}_1.fastq"
  renamed_fastq2="$FASTQ_DIR/${sample_name}_${run_name}_2.fastq"
  renamed_single_fastq="$FASTQ_DIR/${sample_name}_${run_name}.fastq"

  [ -f "$fastq1" ] && mv "$fastq1" "$renamed_fastq1"
  [ -f "$fastq2" ] && mv "$fastq2" "$renamed_fastq2"
  [ -f "$single_fastq" ] && mv "$single_fastq" "$renamed_single_fastq"

  echo "Completed downloading and renaming files for $sample_name ($run_name)"
done < "$INPUT_FILE"

