
# DQA1, HLA-DQ
# for each of these, try to grab gwas.api_query()
# if we return anything from gwas.api_query(f"gwasinfo"):
# if trait and sample_size: 
# trait, sample_size

def protein():
    pass

if __name__ == "__main__":
    protein()

import ieugwaspy as gwas
import requests
import os
import json

creds = json.load(open(".ieugwaspy.json"))
token = creds.get("jwt")

headers = {"Authorization", token}

genes = ["DQA1", "HLA-DQ", "IL2RA"]
for gene in genes:
    try:
        result = requests.get(
        f"https://api.opengwas.io/api/gwas?q={gene}",
        headers=headers
    )
        if result and result.get("trait") and result.get("sample_size"):
            print(f"{gene}: {result['trait']}, n={result['sample_size']}")
    except Exception as e:
        print(f"{gene}: query failed — {e}")