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
    # Load any previously summarized papers (for persistent tracking)
    try:
        with open('../data/previous_papers.json') as f:
            previous_papers = json.load(f)
    except FileNotFoundError:
        previous_papers = []

    # Extract known PMIDs to prevent duplicates
    seen_pmids = {p.get("id", "").split("/")[-2] for p in previous_papers if "id" in p}

    # Load new batch of papers from fetch_pubmed
    with open('../data/new_papers.json') as f:
        new_papers = json.load(f)

    # Filter out already-summarized PMIDs
    unseen_papers = [
        p for p in new_papers
        if p.get("id", "").split("/")[-2] not in seen_pmids
    ]

    # Summarize only new ones
    for paper in unseen_papers:
        summary = summarize_paper(paper.get('title', ''), paper.get('abstract', ''))
        paper['summary'] = summary

    # Append the newly summarized papers to the previous list
    updated_papers = previous_papers + unseen_papers

    # Save the combined list back to previous_papers.json
    with open('../data/previous_papers.json', 'w') as f:
        json.dump(updated_papers, f, indent=2)
