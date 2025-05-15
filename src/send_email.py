import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from configparser import ConfigParser
import json

# Load config
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

# Load the most recent summaries from previous_papers.json
try:
    with open('../data/previous_papers.json') as f:
        papers = json.load(f)
except FileNotFoundError:
    print("❌ No previous_papers.json file found.")
    exit(1)

# Filter for those just summarized (assumes they have 'summary' but aren't in email history)
new_summaries = [p for p in papers if 'summary' in p]

if not new_summaries:
    print("🟡 No new summaries found to email.")
    exit(0)

# Build the email body
body = "\n\n".join([
    f"Title: {p.get('title', 'N/A')}\n"
    f"Journal: {p.get('journal', 'Unknown')}\n"
    f"Link: {p.get('id', 'N/A')}\n"
    f"Summary: {p.get('summary', 'N/A')}"
    for p in new_summaries
])

send_email(subject="New Microbiome Literature Summary", body=body)
print("✅ Email sent successfully.")
