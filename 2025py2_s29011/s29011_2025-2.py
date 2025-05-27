#!/usr/bin/env python3

import csv
import sys
import time
from datetime import datetime
from typing import List, Optional

from Bio import Entrez, SeqIO
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


class NCBIRetriever:

    def __init__(self, email: str, api_key: str):
        Entrez.email = email
        Entrez.api_key = api_key

        self.webenv: Optional[str] = None
        self.query_key: Optional[str] = None
        self.count: int = 0

    def search(self, taxid: str, min_len: Optional[int], max_len: Optional[int]) -> int:
        query = f"txid{taxid}[Organism]"
        if min_len is not None or max_len is not None:
            _min = min_len if min_len is not None else 1
            _max = max_len if max_len is not None else 1_000_000_000
            query += f" AND {_min}:{_max}[SLEN]"

        with Entrez.esearch(db="nucleotide", term=query, usehistory="y", retmax=0) as h:
            res = Entrez.read(h)

        self.webenv = res["WebEnv"]
        self.query_key = res["QueryKey"]
        self.count = int(res["Count"])
        return self.count

    def _fetch_batch(self, start: int, size: int = 500):
        if self.webenv is None or self.query_key is None:
            raise RuntimeError("error")

        with Entrez.efetch(
            db="nucleotide",
            rettype="gb",
            retmode="text",
            retstart=start,
            retmax=size,
            webenv=self.webenv,
            query_key=self.query_key,
        ) as h:
            yield from SeqIO.parse(h, "gb")

    def fetch_all(self, limit: Optional[int] = None, delay: float = 0.34):
        remaining = self.count if limit is None else min(limit, self.count)
        start = 0
        while remaining > 0:
            batch = min(remaining, 500)
            for rec in self._fetch_batch(start, batch):
                yield rec
            remaining -= batch
            start += batch
            time.sleep(delay)

def write_csv(records: List[SeqIO.SeqRecord], path: str):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Accession", "Length", "Description"])
        for r in records:
            w.writerow([r.id, len(r.seq), r.description])


def plot_lengths(records: List[SeqIO.SeqRecord], path: str):
    sorted_recs = sorted(records, key=lambda r: len(r.seq), reverse=True)
    accs = [r.id for r in sorted_recs]
    lengths = [len(r.seq) for r in sorted_recs]

    plt.figure(figsize=(max(8, len(accs) * 0.3), 6))
    plt.plot(accs, lengths, marker="o")
    plt.xticks(rotation=90, fontsize=7)
    plt.xlabel("GenBank accession")
    plt.ylabel("Sequence length (bp)")
    plt.title("Sequence lengths sorted descending")
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def prompt_int(msg: str) -> Optional[int]:
    val = input(msg).strip()
    if val == "":
        return None
    try:
        return int(val)
    except ValueError:
        print("Wartość musi być liczbą całkowitą")
        return None


def main():

    email = input("Enter your email address for NCBI: ").strip()
    api_key = input("Enter your NCBI API key: ").strip()
    taxid = input("Enter taxonomic ID (taxid) of the organism: ").strip()

    min_len = prompt_int("Minimal sequence length ")
    max_len = prompt_int("Maximal sequence length ")
    max_records = prompt_int("Maximum records to download ")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    prefix = f"taxid_{taxid}_{timestamp}"
    csv_file = f"{prefix}.csv"
    png_file = f"{prefix}.png"

    retriever = NCBIRetriever(email, api_key)

    total = retriever.search(taxid, min_len, max_len)
    if total == 0:
        print("Not found")
        sys.exit(0)

    records = list(retriever.fetch_all(limit=max_records))

    if not records:
        print("Download failed")
        sys.exit(1)

    write_csv(records, csv_file)
    print(f"CSV report saved -> {csv_file}")

    plot_lengths(records, png_file)
    print(f"Length plot saved -> {png_file}\n")

if __name__ == "__main__":
    main()
