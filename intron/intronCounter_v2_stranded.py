#!/usr/bin/python3
"""Intron counter."""
import argparse

import pybedtools
import pysam


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Count number of unprocessed reads in each intron.")
    parser.add_argument("-b", "--bam", help="Input BAM file.", required=True)
    parser.add_argument("-g", "--gtf", help="Input GTF file with introns.", required=True)
    parser.add_argument("-o", "--output", help="Output file.", required=True)
    return parser.parse_args()


def is_read_on_correct_strand(strand, read):
    """Check if read is on expected strand."""
    if strand == "+" and ((read.is_read1 and read.is_reverse) or (read.is_read2 and not read.is_reverse)):
        return True
    elif strand == "-" and ((read.is_read1 and not read.is_reverse) or (read.is_read2 and read.is_reverse)):
        return True

    return False


def intronCounter(file, gtf_file, output_file):
    """Count number of reads in intron."""
    with pysam.AlignmentFile(file, "rb") as bamfile, open(output_file, "wt") as fileout:
        for intron in pybedtools.BedTool(gtf_file):
            fileout.write(str(intron).rstrip("\n") + "\t")

            if len(intron) <= 10:
                fileout.write("0" + "\n")
                continue

            intron_count = 0
            span_count = 0
            for read in bamfile.fetch(intron.chrom, intron.start + 5, intron.end - 5):
                if read.is_duplicate:
                    continue
                if read.is_qcfail:
                    continue
                if read.is_secondary:
                    continue
                if not read.is_proper_pair:
                    continue

                # This check is for ISR reads:
                # Filtering by correct strand is done for ISR type (Illumina truseq):
                # https://salmon.readthedocs.io/en/latest/library_type.html
                # TODO: Generalize this for different types of kits.
                if is_read_on_correct_strand(intron.strand, read):

                    span_count += 1
                    if "N" not in read.cigarstring:
                        intron_count += 1

            # Calculate rpkm
            rpkm = (float(intron_count) * 10 ** 9) / (bamfile.mapped * len(intron))

            fileout.write("\t".join(map(str, [intron_count, span_count, rpkm])) + "\n")


if __name__ == "__main__":
    args = parse_arguments()
    intronCounter(args.bam, args.gtf, args.output)
