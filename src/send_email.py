import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from configparser import ConfigParser
import json
import sys

# Load configuration
config = ConfigParser()
config.read('../config/config.ini')

def send_email(subject, body):
    msg = MIMEMultipart()
    msg['From'] = config['email']['from_email']
    msg['To'] = config['email']['to_email']
    msg['Subject'] = subject
    msg.attach(MIMEText(body, 'plain'))

    server = smtplib.SMTP(config['email']['smtp_server'], int(config['email']['smtp_port']))
    server.starttls()
    server.login(config['email']['username'], config['email']['password'])
    server.sendmail(msg['From'], msg['To'], msg.as_string())
    server.quit()

# Load previously summarized papers
try:
    with open('../data/previous_papers.json') as f:
        all_papers = json.load(f)
    print(f"✅ Loaded {len(all_papers)} previous summarized papers.")
except FileNotFoundError:
    print("❌ previous_papers.json not found.")
    sys.exit(1)

# Load the list of just-summarized IDs
try:
    with open('../data/just_summarized_ids.json') as f:
        just_ids = set(json.load(f))
    print(f"✅ Loaded {len(just_ids)} just-summarized paper IDs.")
except FileNotFoundError:
    print("❌ just_summarized_ids.json not found.")
    sys.exit(1)

if not just_ids:
    print("🟡 No new papers to email. Skipping.")
    sys.exit(0)

# Filter to only the current run’s summaries
to_send = [p for p in all_papers if p.get("id") in just_ids]

if not to_send:
    print("🟡 No matching summaries found to email. Skipping.")
    sys.exit(0)

# Build the email content
body = "\n\n".join([
    f"Title: {p.get('title', 'N/A')}\n"
    f"Journal: {p.get('journal', 'Unknown')}\n"
    f"Link: {p.get('id', 'N/A')}\n"
    f"Summary: {p.get('summary', 'N/A')}"
    for p in to_send
])

send_email(subject="New Microbiome Literature Summary", body=body)
print(f"✅ Email sent with {len(to_send)} papers.")
