#!/usr/bin/bash
# Script to monitor and resubmit Kallisto mapping jobs
KALLISTO_SCRIPT="kallisto_mapping_sortOlder.sh"
INTERVAL=1600  # Default interval (seconds)
LOG_FILE="kallisto_monitor.log"

# Function to handle cleanup on exit
cleanup() {
  echo "$(date): Exiting Kallisto monitoring script." >> $LOG_FILE
  exit 0
}
trap cleanup SIGINT SIGTERM

# Infinite loop to check and resubmit the job
while true; do
  echo "$(date): Checking for active kallisto_mapping jobs..." | tee -a $LOG_FILE

  # Check if any kallisto_mapping job is present in the queue (any state)
  running_jobs=$(squeue -u $USER --name=kallisto_mapping | grep -v JOBID | wc -l)
  
  if [ "$running_jobs" -gt 0 ]; then
    echo "$(date): Kallisto mapping job is still in the queue. Current jobs:" | tee -a $LOG_FILE
    squeue -u $USER --name=kallisto | tee -a $LOG_FILE
    echo "$(date): Waiting 10 minutes before re-checking..." | tee -a $LOG_FILE
    sleep 600
    continue
  fi
  
  echo "$(date): No active or queued kallisto_mapping jobs. Submitting a new job..." | tee -a $LOG_FILE
  sbatch $KALLISTO_SCRIPT
  
  echo "$(date): Waiting for $INTERVAL seconds before the next check..." | tee -a $LOG_FILE
  sleep $INTERVAL
done

