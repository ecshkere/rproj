import requests
import pandas as pd
import re
import time
from bs4 import BeautifulSoup
import sys

input_path = sys.argv[1]
output_path = sys.argv[2] if len(sys.argv) > 2 else "uniprot_output.csv"

def fetch_variants(uniprot_id):
    url = f"https://www.ebi.ac.uk/proteins/api/variation/{uniprot_id}"
    headers = {"Accept": "application/json"}
    r = requests.get(url, headers=headers)
    r.raise_for_status()
    data = r.json()

    rows = []
    for feature in data.get("features", []):
        pos = feature.get("begin")
        for xref in feature.get("xrefs", []):
            if xref.get("name") in {"dbSNP", "ExAC", "ESP", "gnomAD", "TOPMed", "1000Genomes"}:
                rsid = xref.get("id")
                rows.append({
                    "uniprot_id": uniprot_id,
                    "rsid": rsid,
                    "protein_pos": pos
                })
    return rows

def fetch_domains(uniprot_id):
    url = f"https://web.expasy.org/cgi-bin/protparam/protparam?{uniprot_id}"
    r = requests.get(url)
    if r.status_code != 200:
        print(f"Failed {uniprot_id}, status {r.status_code}")
        return []

    soup = BeautifulSoup(r.text, "html.parser")
    pre = soup.find("pre")
    if not pre:
        return []

    domains = []
    for line in pre.get_text().splitlines():
        # Example line: "FT   DOMAIN         74-145  KRAB"
        m = re.match(r"FT\s+(\w+)\s+(\d+)-(\d+)\s+(.+)", line.strip())
        if m:
            feature, start, end, name = m.groups()
            if feature in ["CHAIN", "COMPBIAS"]:
                continue
            domains.append({
                "feature": feature,
                "start": int(start),
                "end": int(end),
                "name": name.strip()
            })
    return domains

def map_position_to_domain(position, domains):
    matches = []
    for d in domains:
        if d["start"] <= position <= d["end"]:
            matches.append(f"{d['feature']}|{d['name']}")
    return "; ".join(matches) if matches else "NA"

df = pd.read_csv(input_path, sep="\t")
id_map = df.groupby("Entry")["rsid"].apply(list).to_dict()

all_rows = []
for uid, rsids in id_map.items():
    rows = fetch_variants(uid)
    filtered = [r for r in rows if r["rsid"] in rsids]
    all_rows.extend(filtered)

df = pd.DataFrame(all_rows)
df = df.drop_duplicates()
df.reset_index(drop=True, inplace=True)

domain_cache = {}  # cache results to avoid duplicate requests

results = []
for i, row in df.iterrows():
    uniprot_id = row["uniprot_id"]
    pos = int(row["protein_pos"])

    if uniprot_id not in domain_cache:
        domain_cache[uniprot_id] = fetch_domains(uniprot_id)
        time.sleep(0.5)  # polite delay

    domains = domain_cache[uniprot_id]
    mapped = map_position_to_domain(pos, domains)

    results.append(mapped)
    print(uniprot_id, "done", i + 1, "/", df.shape[0]) # FIX

df["Domain"] = results

# Save new table
df.sort_values(by='Domain', inplace=True)
df.to_csv(output_path, sep="\t", index=False)
