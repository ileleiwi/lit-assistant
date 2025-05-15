import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from configparser import ConfigParser
import json

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

if __name__ == "__main__":
    with open('../data/previous_papers.json') as f:
        summaries = json.load(f)

    body = "\n\n".join([
    f"Title: {p.get('title', 'N/A')}\n"
    f"Journal: {p.get('journal', 'Unknown')}\n"
    f"Link: {p.get('id', 'N/A')}\n"
    f"Summary: {p.get('summary', 'N/A')}"
    for p in summaries
    ])



    send_email(subject="Daily Microbiome Literature Update", body=body)

