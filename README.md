# Literature Assistant

Automated AI-powered assistant to scan microbiome literature, summarize new papers, and deliver updates via email.

## Project Structure

```
lit-assistant/
├── config/ # Configurations (ignored by git)
├── data/ # Paper data & summaries (ignored by git)
├── logs/ # Logs for debugging (ignored by git)
├── src/ # Python scripts for workflow
│ ├── fetch_papers.py
│ ├── summarize.py
│ ├── send_email.py
│ └── daily_lit_summary.py
├── environment.yml # Conda environment file
└── run_daily.sh # Cron script for daily execution
```


## Setup Instructions
- Clone repository
- Create conda environment from `environment.yml`
- Populate `config/config.ini`
- Run `run_daily.sh` via cron job

## Usage
Automated via cron. Outputs emailed daily.

## Security
Never commit your config files (`config.ini`) or sensitive credentials.



