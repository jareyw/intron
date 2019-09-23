"""Intron counter."""
import argparse
import csv
import gzip

import pybedtools
import pysam


def parse_args():
    """Parse arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bam", help="Input bam file.")
    parser.add_argument("gtf", help="Custom GTF file (should only contain introns features).")
    parser.add_argument("outfile", help="Name of output file (should have .tsv.gz extension).")
    parser.add_argument("-o", "--min-overhang", help="Minimal overhang of read over intron.", default=5)

    return parser.parse_args()


def calculate_rpkm(reads, feature_size, mapped_reads):
    """Compute RPKM - Reads Per Kb of feature, per Million mapped reads."""
    return reads * 10 ** 9 / (mapped_reads * feature_size)


def count_introns(bam, gtf, outfile, min_overhang=5):
    """Count introns for single bam file."""
    if not outfile.endswith(".tsv.gz"):
        raise ValueError("Output file should have an .tsv.gz extension.")

    with gzip.open(outfile, "wt") as out, pysam.AlignmentFile(bam, mode="rb", require_index=True) as bamfile:

        csvwriter = csv.writer(out, delimiter="\t", lineterminator="\n")

        # The GTF file only contains introns so it makes sense to refer
        # to each segment in GTF as "intron".
        for intron in pybedtools.BedTool(gtf):

            start = intron.start + min_overhang
            end = intron.end - min_overhang
            if start >= end:
                csvwriter.writerow(intron[:] + ["0", "0", "0.0"])
                continue

            intron_count = 0
            span_count = 0
            for read in bamfile.fetch(intron.chrom, start, end):

                if read.is_qcfail:
                    continue
                if read.is_duplicate:
                    continue
                if read.is_secondary:
                    continue
                if not read.is_proper_pair:
                    continue
                read_strand = "-" if read.is_reverse else "+"
                if read_strand != intron.strand:
                    continue

                if "N" not in read.cigarstring:
                    intron_count += 1
                else:
                    span_count += 1

            # Calculate rpkm per genomic locus for every intron
            rpkm = calculate_rpkm(intron_count, len(intron), bamfile.mapped)

            csvwriter.writerow(intron[:] + [intron_count, span_count, rpkm])


def main():
    """Invoke when run directly as a program."""
    args = parse_args()
    count_introns(bam=args.bam, gtf=args.gtf, outfile=args.outfile, min_overhang=args.min_overhang)


if __name__ == "__main__":
    main()
