import openai
from configparser import ConfigParser

config = ConfigParser()
config.read('../config/config.ini')
openai.api_key = config['openai']['api_key']
model = config['openai']['model']

def summarize_paper(title, abstract):
    prompt = (
        "You are summarizing microbiome research papers succinctly. "
        "Highlight the novelty, methods, and potential implications.\n\n"
        f"Title: {title}\nAbstract: {abstract}\n\nSummary:"
    )

    response = openai.ChatCompletion.create(
        model=model,
        messages=[{"role": "user", "content": prompt}],
        temperature=0.2,
        max_tokens=200,
    )
    return response.choices[0].message.content.strip()

if __name__ == "__main__":
    import json
    with open('../data/previous_papers.json') as f:
        papers = json.load(f)

    summaries = []
    for paper in papers:
        summary = summarize_paper(paper['title'], paper['abstract'])
        summaries.append({
            'title': paper['title'],
            'link': paper['id'],
            'summary': summary
        })

    with open('../data/previous_papers.json', 'w') as f:
        json.dump(summaries, f, indent=2)

