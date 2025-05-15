import json
import openai
from configparser import ConfigParser
from datetime import datetime

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
#     response = openai.ChatCompletion.create(...)
#     return response.choices[0].message.content.strip()

def summarize_paper(title, abstract):
    return f"(MOCK SUMMARY)\nTitle: {title[:60]}...\nAbstract: {abstract[:100]}..."

if __name__ == "__main__":
    try:
        with open('../data/previous_papers.json') as f:
            previous = json.load(f)
    except FileNotFoundError:
        previous = []

    seen_ids = {p.get("id", "") for p in previous if "id" in p}

    try:
        with open('../data/new_papers.json') as f:
            pubmed = json.load(f)
    except FileNotFoundError:
        pubmed = []

    try:
        with open('../data/biorxiv_papers.json') as f:
            biorxiv = json.load(f)
    except FileNotFoundError:
        biorxiv = []

    all_new = pubmed + biorxiv
    new_to_summarize = []
    now = datetime.now().isoformat()

    for p in all_new:
        pid = p.get("id", "")
        if pid not in seen_ids:
            p['summary'] = summarize_paper(p.get("title", ""), p.get("abstract", ""))
            p['summarized_on'] = now
            new_to_summarize.append(p)

    updated = previous + new_to_summarize

    with open('../data/previous_papers.json', 'w') as f:
        json.dump(updated, f, indent=2)

        # Write even if it's an empty list
    just_ids = [p["id"] for p in new_to_summarize if "id" in p]
    with open('../data/just_summarized_ids.json', 'w') as f:
        json.dump(just_ids, f)

    print(f"\n✅ Summarized {len(new_to_summarize)} new papers.")

