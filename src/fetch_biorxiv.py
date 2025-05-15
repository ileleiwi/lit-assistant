import feedparser
import json
import re
import unicodedata
from bs4 import BeautifulSoup

# Expanded and robust keyword list
KEYWORDS = [
    "microbiome", "metagenomics", "qsip", "sip", "16s",
    "multi omics", "multi-omics", "transcriptomics", "metabolomics",
    "proteomics", "soil", "machine learning", "machine-learning",
    "artificial intelligence", "deep learning", "electron microscopy"
]

FEED_URL = 'https://www.biorxiv.org/rss/latest.xml'

def clean_text(text):
    # Remove HTML, normalize characters, standardize hyphens and whitespace
    text = BeautifulSoup(text, "html.parser").get_text()
    text = unicodedata.normalize("NFKD", text)
    text = text.replace("-", " ")
    text = text.lower()
    return re.sub(r'\s+', ' ', text).strip()

def fetch_biorxiv():
    feed = feedparser.parse(FEED_URL)
    new_papers = []

    for entry in feed.entries:
        title = entry.title
        summary = entry.summary
        link = entry.link

        # Prepare cleaned text for keyword matching
        full_text = clean_text(title + " " + summary)

        matched_keywords = [kw for kw in KEYWORDS if kw.lower() in full_text]

        if matched_keywords:
            print(f"✅ Matched bioRxiv paper: {title}")
            print(f"   Keywords matched: {', '.join(matched_keywords)}\n")

            new_papers.append({
                "id": link,
                "title": title,
                "abstract": summary,
                "journal": "bioRxiv"
            })

    return new_papers

if __name__ == "__main__":
    papers = fetch_biorxiv()
    with open('../data/biorxiv_papers.json', 'w') as f:
        json.dump(papers, f, indent=2)

    print(f"🧠 Found {len(papers)} bioRxiv papers matching keywords.")
