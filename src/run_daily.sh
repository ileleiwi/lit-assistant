#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate lit-assistant
cd /path/to/lit-assistant/src
python daily_lit_summary.py

