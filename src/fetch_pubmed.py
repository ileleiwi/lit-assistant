from Bio import Entrez
from configparser import ConfigParser
import xml.etree.ElementTree as ET
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
        xml_data = fetch_handle.read()
        root = ET.fromstring(xml_data)

        for article in root.findall(".//PubmedArticle"):
            try:
                pmid = article.findtext(".//PMID")
                title = article.findtext(".//ArticleTitle", default="No title available")
                abstract = article.findtext(".//AbstractText", default="No abstract available")
                journal = article.findtext(".//Journal/Title", default="Unknown")

                papers.append({
                    "id": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "N/A",
                    "title": title.strip(),
                    "abstract": abstract.strip(),
                    "journal": journal.strip()
                })

            except Exception as e:
                print(f"Error processing article: {e}")
                continue

        fetch_handle.close()

    return papers

if __name__ == "__main__":
    papers = fetch_pubmed(query, max_results)
    with open('../data/previous_papers.json', 'w') as f:
        json.dump(papers, f, indent=2)
