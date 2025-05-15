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
        print(f"❌ {name} failed. See logs.")
    else:
        logging.info(f"{name} completed successfully.")
        print(f"✅ {name} ran successfully.")

def count_papers(path):
    try:
        with open(path) as f:
            papers = json.load(f)
        return len(papers)
    except:
        return 0

if __name__ == "__main__":
    print("🚀 Starting literature assistant pipeline...\n")

    run_script('fetch_pubmed.py')
    run_script('fetch_biorxiv.py')
    run_script('summarize.py')

    # Count how many new papers exist before emailing
    num_pubmed = count_papers('../data/new_papers.json')
    num_biorxiv = count_papers('../data/biorxiv_papers.json')
    total_new = num_pubmed + num_biorxiv

    print(f"\n📊 Found {num_pubmed} new PubMed papers and {num_biorxiv} new bioRxiv papers.")

    if total_new > 0:
        run_script('send_email.py')
    else:
        print("🟡 No new papers to email. Skipping send_email.py.")
        logging.info("No new papers found — skipping send_email.py.")
