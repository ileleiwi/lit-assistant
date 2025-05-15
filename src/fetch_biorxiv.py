import feedparser
import json

# Define keywords you care about
KEYWORDS = ['microbiome', 'metagenomics', 'qsip', 'multi-omics', 'machine learning', 'soil']

# RSS feed for bioRxiv latest preprints
FEED_URL = 'https://www.biorxiv.org/rss/latest.xml'

def fetch_biorxiv():
    feed = feedparser.parse(FEED_URL)
    new_papers = []

    for entry in feed.entries:
        title = entry.title
        summary = entry.summary
        link = entry.link

        # Simple keyword match (case-insensitive)
        if any(kw.lower() in (title + summary).lower() for kw in KEYWORDS):
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

    print(f"Found {len(papers)} bioRxiv papers matching keywords.")
