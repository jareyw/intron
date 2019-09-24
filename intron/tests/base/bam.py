"""Test utilities for GTF files."""

import os
import random

import pysam

from .fasta import FastaTestCaseMixin
from .fastq import FastqTestCaseMixin

# Mapping between SAM bits and attributes of pysam.AlignedSegment
# This mapping originates from here:
# https://github.com/pysam-developers/pysam/blob/master/htslib/htslib/sam.h#L114
ATTRIBUTE_BIT_MAP = {
    # 1 0x1 template having multiple segments in sequencing
    "is_paired": 1,
    # 2 0x2 each segment properly aligned according to the aligner
    "is_proper_pair": 2,
    # 4 0x4 segment unmapped
    "is_unmapped": 4,
    # 8 0x8 next segment in the template unmapped
    "mate_is_unmapped": 8,
    # 16 0x10 SEQ being reverse complemented
    "is_reverse": 16,
    # 32 0x20 SEQ of the next segment in the template being reverse complemented
    "mate_is_reverse": 32,
    # 64 0x40 the first segment in the template
    "is_read1": 64,
    # 128 0x80 the last segment in the template
    "is_read2": 128,
    # 256 0x100 secondary alignment
    "is_secondary": 256,
    # 512 0x200 not passing filters, such as platform/vendor quality controls
    "is_qcfail": 512,
    # 1024 0x400 PCR or optical duplicate
    "is_duplicate": 1024,
    # 2048 0x800 supplementary alignment
    "is_supplementary": 2048,
}


class BamTestCaseMixin:
    """Mixin for manipulation of BAM files in tests."""

    def get_flag_value(self, **kwargs):
        """Assign bitwise flag."""
        flag = 0
        for attribute, bit in ATTRIBUTE_BIT_MAP.items():
            if kwargs.get(attribute, False):
                flag += bit

        return flag

    def make_bam_segment(
        self,
        qname=None,
        flag=0,
        rname=0,
        pos=0,
        mapq=20,
        cigar=None,
        rnext=0,
        pnext=0,
        tlen=0,
        seq=None,
        qual=None,
        tags=None,
        **kwargs,
    ):
        """
        Return pysam.AlignedSegment object.

        Each pysam.AlignedSegment element has 11 mandatory tab-separated
        fields. Additionaly there is TAGS field for additional info::

            QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL TAGS

        Columns RNEXT and RNAME are left undefined, since we are only interested
        in sigle-end reads. TLEN is also set to undefined (value 0).

        SEQ and QUAL are randomly generated, their size is determined by cigar value.

        Read more about SAM file specifications here:
        [1] https://samtools.github.io/hts-specs/SAMv1.pdf

        """
        segment = pysam.AlignedSegment()

        if qname is None:
            qname = "read-{}".format(random.randrange(1000, 9999))
        segment.query_name = qname

        segment.flag = flag or self.get_flag_value(**kwargs)
        segment.reference_id = rname
        segment.reference_start = pos
        segment.mapping_quality = mapq
        segment.cigar = cigar

        segment.next_reference_id = rnext
        segment.next_reference_start = pnext
        segment.template_length = tlen

        length = sum([length for (operation, length) in segment.cigartuples if operation in [0, 1, 4, 7, 8]])
        if seq is None:
            seq = FastaTestCaseMixin.make_fasta_sequence(size=length, include_n=False)
        segment.query_sequence = seq

        if qual is None:
            qual = pysam.qualitystring_to_array(FastqTestCaseMixin.make_fastq_sequence(size=length))
        segment.query_qualities = qual

        if tags is not None:
            segment.tags = tags.items()

        return segment

    def make_bam(self, chroms, segments):
        """
        Make a BAM file.

        Each BAM file consists of header and alignment entries. Main
        content of the header are the reference sequence definitions
        (chromosomes). Alignment entries (segments) are records of how
        each read is mapped to the reference.

        The structure of parameter `chroms` should be the following::

            [('chr1', 3000), ('chr2', 2000)],

        Parameter `segments` should be a list of pysam.AlignedSedment objects.

        Links for further read:
        https://samtools.github.io/hts-specs/SAMv1.pdf
        http://pysam.readthedocs.io/en/stable/usage.html#creating-bam-cram-sam-files-from-scratch

        """
        fname = self.get_filename(extension="bam")

        chromosomes = [{"SN": chrom, "LN": size} for chrom, size in chroms]
        header = {"HD": {"VN": "1.0"}, "SQ": chromosomes}

        with pysam.AlignmentFile(fname, "wb", header=header) as outf:
            for segment in segments:
                if isinstance(segment, dict):
                    segment = self.make_bam_segment(**segment)
                outf.write(segment)

        pysam.index(fname)

        return os.path.abspath(fname)
