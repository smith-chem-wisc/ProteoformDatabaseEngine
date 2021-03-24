# see https://www.ebi.ac.uk/training/online/sites/ebi.ac.uk.training.online/files/UniProt_programmatically_py3.pdf

import requests
import sys
import yaml

# find proteome
BASE = 'http://www.uniprot.org'
KB_ENDPOINT = '/proteomes/'
TOOL_ENDPOINT = '/uploadlists/'

with open("config/config.yaml", 'r') as stream:
   data = yaml.safe_load(stream)

query = data["species"] # read config
organism = data["organism"].lower()

# special case
if query == 'canis_familiaris':
    query = 'canis_lupus_familiaris'

payload = {
    'query': query,
    'sort': 'score',
    'format': 'tab',
    }

proteome_res = requests.get(BASE + KB_ENDPOINT, params=payload, stream=True)
proteome_res.raise_for_status() # throw an error for bad status code

results = proteome_res.text.split('\n')[1:]

for r in results:
    splt = r.split('\t')
    if organism in splt[1].lower():
        proteome = splt[0]
        break
