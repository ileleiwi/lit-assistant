import json
import openai
from configparser import ConfigParser

config = ConfigParser()
config.read('../config/config.ini')
openai.api_key = config['openai']['api_key']
model = config['openai']['model']

# def summarize_paper(title, abstract):
#     prompt = (
#         "You are summarizing microbiome research papers succinctly. "
#         "Highlight the novelty, methods, and potential implications.\n\n"
#         f"Title: {title}\nAbstract: {abstract}\n\nSummary:"
#     )
# 
#     response = openai.ChatCompletion.create(
#         model=model,
#         messages=[{"role": "user", "content": prompt}],
#         temperature=0.2,
#         max_tokens=200,
#     )
#     return response.choices[0].message.content.strip()

def summarize_paper(title, abstract):
    return f"(MOCK SUMMARY)\nTitle: {title[:60]}...\nAbstract: {abstract[:100]}..."

if __name__ == "__main__":
    # Load previously summarized papers
    try:
        with open('../data/previous_papers.json') as f:
            previous_papers = json.load(f)
    except FileNotFoundError:
        previous_papers = []

    seen_pmids = {p.get("id", "").split("/")[-2] for p in previous_papers if "id" in p}

    # Load new PubMed papers
    try:
        with open('../data/new_papers.json') as f:
            pubmed_papers = json.load(f)
    except FileNotFoundError:
        pubmed_papers = []

    # Load bioRxiv papers (optional)
    try:
        with open('../data/biorxiv_papers.json') as f:
            biorxiv_papers = json.load(f)
    except FileNotFoundError:
        biorxiv_papers = []

    # Combine all sources
    combined = pubmed_papers + biorxiv_papers

    # Filter out already summarized papers by ID
    new_to_summarize = [
        p for p in combined
        if p.get("id", "").split("/")[-2] not in seen_pmids
    ]

    # Summarize new papers
    for paper in new_to_summarize:
        summary = summarize_paper(paper.get('title', ''), paper.get('abstract', ''))
        paper['summary'] = summary

    # Append new papers to previous and save
    updated_papers = previous_papers + new_to_summarize

    with open('../data/previous_papers.json', 'w') as f:
        json.dump(updated_papers, f, indent=2)
