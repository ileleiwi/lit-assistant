import requests
import json

KEYWORDS = [
    "microbiome", "metagenomics", "qsip", "sip", "16s",
    "multi omics", "multi-omics", "transcriptomics", "metabolomics",
    "proteomics", "soil", "machine learning", "machine-learning",
    "artificial intelligence", "deep learning", "electron microscopy"
]

def fetch_biorxiv_from_epmc():
    query = ' OR '.join([f'"{kw}"' for kw in KEYWORDS])
    full_query = f'({query}) AND SRC:"biorxiv"'

    url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search'
    params = {
        'query': full_query,
        'format': 'json',
        'pageSize': 25
    }

    response = requests.get(url, params=params)
    data = response.json()

    papers = []
    for hit in data.get("resultList", {}).get("result", []):
        papers.append({
            "id": hit.get("id", ""),
            "title": hit.get("title", ""),
            "abstract": hit.get("abstractText", ""),
            "journal": "bioRxiv",
            "source": "EuropePMC"
        })

    return papers

if __name__ == "__main__":
    papers = fetch_biorxiv_from_epmc()
    with open('../data/biorxiv_papers.json', 'w') as f:
        json.dump(papers, f, indent=2)

    print(f"✅ Found {len(papers)} bioRxiv papers matching keywords from Europe PMC.")
