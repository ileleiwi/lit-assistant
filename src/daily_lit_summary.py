import subprocess
import os
import logging
import json

logging.basicConfig(filename='../logs/lit_assistant.log', level=logging.INFO, format='%(asctime)s %(message)s')

def run_script(name):
    script = os.path.join(os.path.dirname(__file__), name)
    result = subprocess.run(['python', script], capture_output=True, text=True)
    if result.returncode != 0:
        logging.error(f"{name} failed: {result.stderr}")
    else:
        logging.info(f"{name} completed successfully.")

if __name__ == "__main__":
    # Step 1: Fetch PubMed and bioRxiv papers
    run_script('fetch_pubmed.py')
    run_script('fetch_biorxiv.py')

    # Step 2: Summarize papers (merges and updates previous_papers.json)
    run_script('summarize.py')

    # Step 3: Only send email if there are new papers
    try:
        with open('../data/new_papers.json') as f:
            new_pubmed = json.load(f)
    except FileNotFoundError:
        new_pubmed = []

    try:
        with open('../data/biorxiv_papers.json') as f:
            new_biorxiv = json.load(f)
    except FileNotFoundError:
        new_biorxiv = []

    total_new = len(new_pubmed) + len(new_biorxiv)

    if total_new > 0:
        run_script('send_email.py')
    else:
        print("🟡 No new papers to email. Skipping send_email.py.")
        logging.info("No new papers found — skipping send_email.py.")
