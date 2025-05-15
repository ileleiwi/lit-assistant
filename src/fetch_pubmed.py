from Bio import Entrez
from configparser import ConfigParser
import xml.etree.ElementTree as ET

config = ConfigParser()
config.read('../config/config.ini')

Entrez.email = config['search']['email']
query = config['search']['query']
max_results = int(config['search']['max_results'])

def fetch_pubmed(query, max_results=10):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="pub+date")
    record = Entrez.read(handle)
    handle.close()
    ids = record["IdList"]

    papers = []
    if ids:
        fetch_handle = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="xml")
        articles = Entrez.read(fetch_handle)
        for article in articles["PubmedArticle"]:
            article_data = article["MedlineCitation"]["Article"]
            title = article_data["ArticleTitle"]
            abstract = article_data.get("Abstract", {}).get("AbstractText", [""])[0]
            pmid = article["MedlineCitation"]["PMID"]
            papers.append({
                "id": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                "title": title,
                "abstract": abstract
            })
        fetch_handle.close()

    return papers

if __name__ == "__main__":
    papers = fetch_pubmed(query, max_results)
    import json
    with open('../data/previous_papers.json', 'w') as f:
        json.dump(papers, f, indent=2)
