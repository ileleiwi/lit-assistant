import requests
import json
from datetime import date, timedelta

KEYWORDS = [
    "microbiome", "metagenomics", "qsip", "sip", "16s",
    "multi omics", "multi-omics", "transcriptomics", "metabolomics",
    "proteomics", "soil", "machine learning", "machine-learning",
    "artificial intelligence", "deep learning", "electron microscopy"
]

def fetch_biorxiv():
    end = date.today()
    start = end - timedelta(days=3)
    url = f"https://api.biorxiv.org/details/biorxiv/{start}/{end}/0"
    response = requests.get(url)
    response.raise_for_status()
    results = response.json()

    try:
        with open('../data/previous_papers.json') as f:
            previous = json.load(f)
            seen_ids = {p.get("id", "") for p in previous if "id" in p}
    except FileNotFoundError:
        seen_ids = set()

    papers = []
    for entry in results.get("collection", []):
        combined = f"{entry.get('title', '').lower()} {entry.get('abstract', '').lower()}"
        if any(k in combined for k in KEYWORDS):
            pid = f"https://www.biorxiv.org/content/{entry['doi']}"
            if pid in seen_ids:
                continue
            papers.append({
                "id": pid,
                "title": entry.get("title", ""),
                "abstract": entry.get("abstract", ""),
                "journal": "bioRxiv",
                "date": entry.get("date", ""),
                "source": "bioRxiv"
            })

    return papers

if __name__ == "__main__":
    papers = fetch_biorxiv()
    with open('../data/biorxiv_papers.json', 'w') as f:
        json.dump(papers, f, indent=2)
