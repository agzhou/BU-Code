#!/bin/bash
#
# Copyright (C) 2001-2025 Verasonics, Inc.
# All worldwide rights and remedies under all intellectual property laws and industrial property laws are reserved.
#
#

# Configurable parameters
TARGET_DIR="/media/verasonics/WD0"   # Change to your NVMe mount point
OUTPUT_CSV="nvme_write_benchmark.csv"
FILE_SIZE_GB=2
BLOCK_SIZE="1M"
RUNTIME_LIMIT=360  # Set to non-zero seconds to auto-stop after this time
INTERVAL=0.3       # Seconds to wait between writes

echo "Starting $RUNTIME_LIMIT second Linux FIO test on drive: $TARGET_DIR for file size: $FILE_SIZE_GB GB"

# Initialize counters
FILE_INDEX=0
TOTAL_WRITTEN_GB=0
START_TIME=$(date +%s)

# Create target directory if it doesn't exist
mkdir -p "$TARGET_DIR"

# Write CSV header
echo "Time_Elapsed(s),Total_Data_Written(GB),Write_Speed(MB/s)" > "$OUTPUT_CSV"

while true; do
    FILE_INDEX=$((FILE_INDEX + 1))
    FILENAME=$(printf "write_%04d.dat" "$FILE_INDEX")
    FILEPATH="$TARGET_DIR/$FILENAME"

    CURRENT_TIME=$(date +%s)
    ELAPSED=$((CURRENT_TIME - START_TIME))

    # Use fio to write a 2GB file with 1M blocks, sequential write
    FIO_OUTPUT=$(fio --name=write_test \
        --filename="$FILEPATH" \
        --filesize="$FILE_SIZE_GB"G \
        --bs=$BLOCK_SIZE \
        --rw=write \
        --ioengine=libaio \
        --direct=1 \
        --verify=0 \
        --numjobs=1 \
        --iodepth=64 \
        --group_reporting=1 \
        --output-format=json)

    # Parse write speed in KB/s from FIO output
    WRITE_SPEED_KB=$(echo "$FIO_OUTPUT" | jq '.jobs[0].write.bw')

    if [[ -z "$WRITE_SPEED_KB" || "$WRITE_SPEED_KB" -eq 0 ]]; then
        echo "Write speed data not available or drive full. Stopping."
        break
    fi

    # Convert to MB/s
    WRITE_SPEED_MB=$(awk "BEGIN {printf \"%.2f\", $WRITE_SPEED_KB / 1024}")

    # Track total data written
    TOTAL_WRITTEN_GB=$((FILE_INDEX*($FILE_SIZE_GB)))

    # Log to CSV
    echo "$ELAPSED,$TOTAL_WRITTEN_GB,$WRITE_SPEED_MB" >> "$OUTPUT_CSV"
    echo "Elapsed time(s): $ELAPSED, Total written(GB): $TOTAL_WRITTEN_GB, Write Speed(MB/s): $WRITE_SPEED_MB"

    # Stop after runtime limit if set
    if [[ $RUNTIME_LIMIT -gt 0 && $ELAPSED -ge $RUNTIME_LIMIT ]]; then
        echo "Reached runtime limit. Exiting."
        break
    fi

    sleep $INTERVAL
done

echo "Benchmark complete. Results saved to $OUTPUT_CSV"
