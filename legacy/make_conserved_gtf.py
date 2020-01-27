#!/usr/bin/python3
"""
Create custom gtf file.

Actions:
- Keep only chr1, ..., chr23, chrX, chrY
- Remove NR noncoding transcripts
- Remove poorly annotated LOC transcripts,
- Keep only exon features
- Remove overlapping genes.
- Remove too long genes.

"""
import argparse
from collections import defaultdict

import pybedtools


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Make a custom GTF file.")
    parser.add_argument("-i", "--input", help="Input GTF file.", required=True)
    parser.add_argument("-o", "--output", help="Output GTF file.", required=True)
    return parser.parse_args()


def process_gtf(input_file, output_file):
    """Parse original GTF."""
    chrom_list = ["chr" + str(i) for i in range(1, 23)]
    chrom_list.extend(["chrX", "chrY"])

    d = defaultdict(lambda: defaultdict(list))
    for segment in pybedtools.BedTool(input_file):
        gene_id = segment.attrs["gene_id"]
        if segment.chrom not in chrom_list:
            continue
        if segment.attrs["transcript_id"].startswith("NR"):
            continue
        if gene_id.startswith("LOC"):
            continue
        if segment[2] != "exon":
            continue

        d[gene_id]["duple"].extend([segment.start, segment.end])
        d[gene_id]["chrom"] = segment.chrom
        d[gene_id]["strand"] = segment.strand

    # For each gene, obtain the broadest coordinates as a [start, stop] "duple"
    gene_ids = list(d.keys())
    for gene_id in gene_ids:
        mincoord = min(d[gene_id]["duple"])
        maxcoord = max(d[gene_id]["duple"])
        if maxcoord - mincoord > 1000000:
            d.pop(gene_id, None)
        d[gene_id]["locus"] = [mincoord, maxcoord]

    # Remove overlapping genes
    overlapping = set()
    for gene in d:
        for gene2 in d:
            if d[gene]["chrom"] != d[gene2]["chrom"]:
                continue
            if d[gene]["strand"] != d[gene2]["strand"]:
                continue
            if gene == gene2:
                continue

            x, y = d[gene]["locus"]
            x2, y2 = d[gene2]["locus"]

            if x2 < y and y2 > x:
                overlapping.add(gene)
                overlapping.add(gene2)
                break
    for gene in overlapping:
        d.pop(gene, None)

    with open(output_file, "w") as handle:
        for segment in pybedtools.BedTool(input_file):
            if segment.attrs["gene_id"] not in d:
                continue
            if not segment.attrs["transcript_id"].startswith("NM"):
                continue
            if segment[2] != "exon":
                continue

            handle.write(str(segment))


if __name__ == "__main__":
    args = parse_arguments()
    process_gtf(args.input, args.output)
