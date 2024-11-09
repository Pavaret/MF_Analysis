#!/bin/bash

# Assign arguments to variables
WORK_DIR=$(pwd)
R_SCRIPT="${WORK_DIR}/Figure5.R"
PYTHON_SCRIPT="${WORK_DIR}/main.py"

# Change to the specified working directory
if [ -d "$WORK_DIR" ]; then
    cd "$WORK_DIR" || exit
    echo "Changed to working directory: $WORK_DIR"
else
    echo "Error: Directory $WORK_DIR does not exist."
    exit 1
fi

# Run the specified R script
if [ -f "$R_SCRIPT" ]; then
    echo "Running R script: $R_SCRIPT"
    Rscript "$R_SCRIPT"
else
    echo "Error: R script $R_SCRIPT not found."
    exit 1
fi

# Run the specified Python script
if [ -f "$PYTHON_SCRIPT" ]; then
    echo "Running Python script: $PYTHON_SCRIPT"
    python "$PYTHON_SCRIPT"
else
    echo "Error: Python script $PYTHON_SCRIPT not found."
    exit 1
fi

echo "Scripts executed successfully."
exit 0
