from Bio import Entrez
from configparser import ConfigParser
import json

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
            try:
                article_data = article["MedlineCitation"]["Article"]
                title = article_data.get("ArticleTitle", "No title available")

                abstract_block = article_data.get("Abstract", {})
                abstract = abstract_block.get("AbstractText", ["No abstract available"])[0]

                pmid = article["MedlineCitation"].get("PMID", "N/A")
                link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid != "N/A" else "N/A"

                journal = article_data.get("Journal", {}).get("Title", "Unknown")

                papers.append({
                    "id": link,
                    "title": title,
                    "abstract": abstract,
                    "journal": journal
                })

            except Exception as e:
                print(f"Skipping article due to error: {e}")
                continue

        fetch_handle.close()

    return papers

if __name__ == "__main__":
    papers = fetch_pubmed(query, max_results)
    with open('../data/previous_papers.json', 'w') as f:
        json.dump(papers, f, indent=2)
