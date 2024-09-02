#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
PYTHON_SCRIPT="$SCRIPT_DIR/delete_finetune.py"

python "$PYTHON_SCRIPT" \
  --config "$PROJECT_ROOT/configs/finetune.yml" \
  --logdir "$PROJECT_ROOT/logs"
