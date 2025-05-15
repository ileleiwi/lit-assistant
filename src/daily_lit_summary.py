import subprocess
import logging
import os

logging.basicConfig(filename='../logs/lit_assistant.log', level=logging.INFO, format='%(asctime)s %(message)s')

def run_script(script_name):
    script_path = os.path.join(os.path.dirname(__file__), script_name)
    result = subprocess.run(['python', script_path], capture_output=True, text=True)
    if result.returncode != 0:
        logging.error(f"{script_name} failed: {result.stderr}")
    else:
        logging.info(f"{script_name} completed successfully.")

if __name__ == "__main__":
    run_script('fetch_pubmed.py')
    run_script('summarize.py')
    run_script('send_email.py')

