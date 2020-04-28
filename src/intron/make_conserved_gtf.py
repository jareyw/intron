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


UCSC hg38: 26485 genes, 52222 transcripts
    NM transcripts: 38793 (133 LOC)
    NR transcripts: 13429 (2085 LOC)

ENSEMBL 92: 60624 genes, 227368 transcripts
    Protein coding transcripts: 83956
    ncRNA: other

    support level 1 transcripts: 29768
    support level 2 transcripts: 37503
    ...

"""
import argparse
from collections import defaultdict

import pybedtools


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Make a custom GTF file.")
    parser.add_argument("-i", "--input", help="Input GTF file.", required=True)
    parser.add_argument("-o", "--output", help="Output GTF file.", required=True)
    parser.add_argument(
        "-s",
        "--source",
        help="Source of GTF file.",
        required=True,
        choices=["UCSC", "ENSEMBL"],
    )
    return parser.parse_args()


def passes_filter_ucsc(segment):
    """Decide whether segment passes filters for UCSC data."""
    # Remove non-coding transcripts
    if not segment.attrs["transcript_id"].startswith("NM"):
        # Manually curated RefSeq transcript identifiers start with NM_ (coding) or NR_ (non-coding)
        # https://genome.ucsc.edu/FAQ/FAQgenes.html
        # Also:
        # https://en.wikipedia.org/wiki/RefSeq
        return False

    if segment.attrs["gene_id"].startswith("LOC"):
        # Remove poorly annotated LOC transcripts
        return False

    return True


def passes_filter_ensembl(segment):
    """Decide whether segment passes filters for ENSEMBL data."""
    # Remove non-coding transcripts
    if not segment.attrs["gene_biotype"] == "protein_coding":
        return False

    support_level = segment.attrs.get("transcript_support_level", "")
    if not support_level.isdigit():
        return False
    if int(support_level) > 2:
        # Remove poorly annotated transcripts
        # https://www.gencodegenes.org/pages/data_format.html
        return False

    return True


def process_gtf(input_file, output_file, source):
    """Parse original GTF."""
    # Ensembl has 1,2,3, ... X, Y chroms.
    chrom_list = [f"{i}" for i in range(1, 23)] + ["X", "Y"]
    if source == "UCSC":
        chrom_list = [f"chr{chr}" for chr in chrom_list]

    d = defaultdict(lambda: defaultdict(list))
    for segment in pybedtools.BedTool(input_file):
        gene_id = segment.attrs["gene_id"]
        if segment.chrom not in chrom_list:
            continue

        if segment[2] != "exon":
            continue

        if source == "UCSC":
            if not passes_filter_ucsc(segment):
                continue
        elif source == "ENSEMBL":
            if not passes_filter_ensembl(segment):
                continue

        d[gene_id]["duple"].extend([segment.start, segment.end])
        d[gene_id]["chrom"] = segment.chrom
        d[gene_id]["strand"] = segment.strand

    # For each gene, obtain the broadest coordinates as a [start, stop] "duple"
    gene_ids = list(d.keys())
    for gene_id in gene_ids:
        mincoord = min(d[gene_id]["duple"])
        maxcoord = max(d[gene_id]["duple"])
        d[gene_id]["locus"] = [mincoord, maxcoord]
        if maxcoord - mincoord > 1000000:
            d.pop(gene_id, None)

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

            if segment[2] != "exon":
                continue

            if source == "UCSC":
                if not passes_filter_ucsc(segment):
                    continue
            elif source == "ENSEMBL":
                if not passes_filter_ensembl(segment):
                    continue

            handle.write(str(segment))


if __name__ == "__main__":
    args = parse_arguments()
    process_gtf(args.input, args.output, args.source)
