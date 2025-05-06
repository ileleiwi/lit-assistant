import arxiv
import json
import os
from configparser import ConfigParser

config = ConfigParser()
config.read('../config/config.ini')

def fetch_arxiv_papers(query, max_results=10):
    search = arxiv.Search(
        query=query,
        max_results=max_results,
        sort_by=arxiv.SortCriterion.SubmittedDate
    )
    return [
        {
            'id': result.entry_id,
            'title': result.title,
            'abstract': result.summary
        } for result in search.results()
    ]

if __name__ == "__main__":
    papers = fetch_arxiv_papers(config['search']['query'], int(config['search']['max_results']))
    with open('../data/previous_papers.json', 'w') as f:
        json.dump(papers, f, indent=2)

