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

if __name__ == "__main__":
    print("🚀 Starting literature assistant pipeline...\n")
    run_script('fetch_pubmed.py')
    run_script('fetch_biorxiv.py')
    run_script('summarize.py')
    run_script('send_email.py')
