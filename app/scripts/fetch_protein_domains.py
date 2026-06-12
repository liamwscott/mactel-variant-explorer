#!/usr/bin/env python3
"""
One-time fetch of protein lengths + Pfam domain annotation for every gene in
data/gene_info.tsv, written to data/protein_domains.tsv.

The Shiny app reads ONLY the bundled TSV at runtime — it never touches the
network. Re-run this script if the gene list changes:

    python3 scripts/fetch_protein_domains.py

Sources (public, no key required):
  * UniProt REST   -> reviewed human accession + canonical sequence length
  * InterPro REST  -> Pfam domains (accession, name, start, end) for that protein

Output columns (one row per domain; genes with no domains get a single
length-only row so the protein backbone can still be drawn):
  Gene_Symbol  UniProt  Protein_Length  Pfam  Domain  Start  End
"""

import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GENE_INFO = os.path.join(ROOT, "data", "gene_info.tsv")
OUT = os.path.join(ROOT, "data", "protein_domains.tsv")

UA = {"User-Agent": "mactel-variant-explorer/1.0 (one-time domain fetch)"}
TIMEOUT = 30
SLEEP = 0.34  # be polite to the public APIs


def get_json(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
        return json.load(r)


def uniprot_lookup(gene):
    """Return (accession, length) for the reviewed human entry, or (None, None)."""
    q = f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true"
    url = ("https://rest.uniprot.org/uniprotkb/search?"
           + urllib.parse.urlencode({
               "query": q,
               "fields": "accession,length",
               "format": "json",
               "size": "1",
           }))
    try:
        d = get_json(url)
    except Exception as e:
        print(f"  ! UniProt error for {gene}: {e}", file=sys.stderr)
        return None, None
    res = d.get("results", [])
    if not res:
        return None, None
    r = res[0]
    return r["primaryAccession"], r["sequence"]["length"]


def pfam_domains(accession):
    """Return list of (pfam_id, name, start, end) for a UniProt accession."""
    url = (f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/"
           f"{accession}?page_size=100")
    try:
        d = get_json(url)
    except urllib.error.HTTPError as e:
        if e.code == 204:  # no content -> no Pfam matches
            return []
        print(f"  ! InterPro error for {accession}: {e}", file=sys.stderr)
        return []
    except Exception as e:
        print(f"  ! InterPro error for {accession}: {e}", file=sys.stderr)
        return []
    out = []
    for r in d.get("results", []):
        m = r["metadata"]
        for prot in r.get("proteins", []):
            for loc in prot.get("entry_protein_locations", []):
                for fr in loc.get("fragments", []):
                    out.append((m["accession"], m["name"],
                                fr["start"], fr["end"]))
    out.sort(key=lambda x: (x[2], x[3]))
    return out


def read_genes(path):
    genes = []
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            g = (row.get("Gene_Symbol") or "").strip()
            if g:
                genes.append(g)
    # de-dup, preserve order
    seen = set()
    return [g for g in genes if not (g in seen or seen.add(g))]


def main():
    genes = read_genes(GENE_INFO)
    print(f"Fetching domains for {len(genes)} genes ...")
    rows = []
    no_uniprot, no_domains = [], []
    for i, gene in enumerate(genes, 1):
        acc, length = uniprot_lookup(gene)
        time.sleep(SLEEP)
        if acc is None:
            no_uniprot.append(gene)
            print(f"[{i:3}/{len(genes)}] {gene:12} -> no UniProt match")
            continue
        doms = pfam_domains(acc)
        time.sleep(SLEEP)
        if not doms:
            no_domains.append(gene)
            rows.append([gene, acc, length, "", "", "", ""])
            print(f"[{i:3}/{len(genes)}] {gene:12} {acc:8} len={length} "
                  f"-> 0 domains")
        else:
            for pf, name, s, e in doms:
                rows.append([gene, acc, length, pf, name, s, e])
            print(f"[{i:3}/{len(genes)}] {gene:12} {acc:8} len={length} "
                  f"-> {len(doms)} domain(s)")

    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Gene_Symbol", "UniProt", "Protein_Length",
                    "Pfam", "Domain", "Start", "End"])
        w.writerows(rows)

    print(f"\nWrote {len(rows)} rows -> {OUT}")
    if no_uniprot:
        print(f"No UniProt match ({len(no_uniprot)}): "
              f"{', '.join(no_uniprot)}")
    if no_domains:
        print(f"No Pfam domains ({len(no_domains)}): "
              f"{', '.join(no_domains)}")


if __name__ == "__main__":
    main()
