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

def count_new_summaries():
    try:
        with open('../data/previous_papers.json') as f:
            all_papers = json.load(f)
        new_papers = [p for p in all_papers if 'summary' in p]
        return len(new_papers)
    except:
        return 0

if __name__ == "__main__":
    print("🚀 Starting literature assistant pipeline...\n")

    run_script('fetch_pubmed.py')
    run_script('fetch_biorxiv.py')
    run_script('summarize.py')

    new_summaries = count_new_summaries()

    print(f"\n📬 Papers summarized in this run: {new_summaries}")

    if new_summaries > 0:
        run_script('send_email.py')
    else:
        print("🟡 No new papers to email. Skipping send_email.py.")
        logging.info("No new papers found — skipping send_email.py.")
