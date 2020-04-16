#!/usr/bin/python3
"""Create custom gtf file."""
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
    """Process GTF file."""
    gene_order = []
    d = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for segment in pybedtools.BedTool(input_file):

        gene_id = segment.attrs["gene_id"]
        transcript_id = segment.attrs["transcript_id"]

        if gene_id not in gene_order:
            gene_order.append(gene_id)

        d[gene_id]["chrom"] = segment.chrom
        d[gene_id]["strand"] = segment.strand
        d[gene_id]["transcript_id"][transcript_id].extend([segment.start, segment.end])

    # sort the lists of exon start/end locations
    for gene_id in d:
        longest = 0
        for transcript_id in d[gene_id]["transcript_id"]:
            d[gene_id]["transcript_id"][transcript_id] = sorted(d[gene_id]["transcript_id"][transcript_id], key=int)
            length = d[gene_id]["transcript_id"][transcript_id][-1] - d[gene_id]["transcript_id"][transcript_id][0]
            if length > longest:
                longest = length
                d[gene_id]["longest_transcript"] = transcript_id
            # Here exons are converted to introns:
            exons = d[gene_id]["transcript_id"][transcript_id]
            introns = []
            for exon1_end, exon2_start in zip(exons[1::2], exons[2::2]):
                introns.append([exon1_end, exon2_start])
            d[gene_id]["transcript_id"][transcript_id] = introns

    f = open(output_file, "w")
    for gene_id in gene_order:
        longest_transcript = d[gene_id]["longest_transcript"]
        # Only keep introns that are shared by all transcripts in gene:
        for duple in d[gene_id]["transcript_id"][longest_transcript]:
            x = [
                "TOSS" if duple not in d[gene_id]["transcript_id"][key] else "KEEP"
                for key in d[gene_id]["transcript_id"].keys()
            ]
            if "TOSS" not in x:
                chrom = d[gene_id]["chrom"]
                start = str(duple[0] + 1)
                end = str(duple[1])
                strand = d[gene_id]["strand"]
                string = "gene_id %s; gene_name %s; transcript_id %s;" % (gene_id, gene_id, longest_transcript)
                f.write("%s\tjarey\tintron\t%s\t%s\t.\t%s\t.\t%s\n" % (chrom, start, end, strand, string))

    f.close()


if __name__ == "__main__":
    args = parse_arguments()
    process_gtf(args.input, args.output)
