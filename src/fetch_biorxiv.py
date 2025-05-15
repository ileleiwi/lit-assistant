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

    new_papers = []

    for entry in results.get("collection", []):
        title = entry.get("title", "").lower()
        abstract = entry.get("abstract", "").lower()
        combined_text = f"{title} {abstract}"
        if any(keyword in combined_text for keyword in KEYWORDS):
            new_papers.append({
                "id": f"https://www.biorxiv.org/content/{entry['doi']}",
                "title": entry.get("title", ""),
                "abstract": entry.get("abstract", ""),
                "journal": "bioRxiv",
                "date": entry.get("date", "")
            })

    return new_papers

if __name__ == "__main__":
    papers = fetch_biorxiv()
    with open("biorxiv_papers.json", "w") as f:
        json.dump(papers, f, indent=2)
    print(f"✅ Saved {len(papers)} matching bioRxiv papers to biorxiv_papers.json")
