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

    # ✅ Load previously emailed PMIDs once
    try:
        with open('../data/previous_papers.json') as f:
            previous = json.load(f)
            seen_pmids = {p.get("id", "").split("/")[-2] for p in previous if "id" in p}
    except FileNotFoundError:
        seen_pmids = set()

    if ids:
        fetch_handle = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="xml")
        xml_data = fetch_handle.read()
        root = ET.fromstring(xml_data)

        for article in root.findall(".//PubmedArticle"):
            try:
                pmid = article.findtext(".//PMID")
                if not pmid or pmid in seen_pmids:
                    continue  # ✅ Skip duplicates or malformed

                title = article.findtext(".//ArticleTitle", default="No title available")
                abstract = article.findtext(".//AbstractText", default="No abstract available")
                journal = article.findtext(".//Journal/Title", default="Unknown")

                print("==== Paper ====")
                print("PMID:", pmid)
                print("Title:", title)
                print("Journal:", journal)
                print("Abstract:", abstract[:60], "...")

                papers.append({
                    "id": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
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

    # ✅ Save new papers to separate file for summarization step
    with open('../data/new_papers.json', 'w') as f:
        json.dump(papers, f, indent=2)
