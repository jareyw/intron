"""Test intron."""
from ngs_test_utils import testcase

from intron import (
    intronCounter_v2_stranded,
    make_conserved_gtf,
    make_conserved_intron_gtf,
)


class IntronTestCase(testcase.NgsTestCase):
    """Test intron."""

    def setUp(self):
        self.outdir = self.get_tmp_dir()
        self.pkw = {"is_paired": True, "is_proper_pair": True}

    def test_intron_ucsc(self):
        gtf = self.make_gtf(
            [
                dict(
                    seqname="chr1",
                    feature="exon",
                    start=99,
                    end=199,
                    strand="+",
                    transcript_id="NM_1",
                    gene_id="G1",
                ),
                dict(
                    seqname="chr1",
                    feature="exon",
                    start=299,
                    end=399,
                    strand="+",
                    transcript_id="NM_1",
                    gene_id="G1",
                ),
            ]
        )
        bam = self.make_bam(
            chroms=[("chr1", 1000)],
            segments=[
                dict(
                    qname="r1",
                    pos=150,
                    pnext=350,
                    cigar=[(0, 75)],
                    is_read1=True,
                    is_reverse=True,
                    **self.pkw,
                ),
                dict(
                    qname="r1",
                    pos=350,
                    pnext=150,
                    cigar=[(0, 75)],
                    is_read2=True,
                    mate_is_reverse=True,
                    **self.pkw,
                ),
            ],
        )

        ann1 = self.get_filename(extension="gtf")
        ann2 = self.get_filename(extension="gtf")
        result = self.get_filename(extension="gtf")

        make_conserved_gtf.process_gtf(gtf, ann1, source="UCSC")
        make_conserved_intron_gtf.process_gtf(ann1, ann2)
        intronCounter_v2_stranded.intronCounter(bam, ann2, result)

        self.assertEqual(
            self.tsv_to_list(result, columns=[3, 4, 9, 10, 11]),
            [["200", "299", "1", "1", "5000000.0"]],
        )

    def test_intron_ensembl(self):
        gtf = self.make_gtf(
            [
                dict(
                    seqname="1",
                    feature="exon",
                    start=99,
                    end=199,
                    strand="+",
                    transcript_id="T1",
                    gene_id="G1",
                    gene_biotype="protein_coding",
                    transcript_support_level="1",
                ),
                dict(
                    seqname="1",
                    feature="exon",
                    start=299,
                    end=399,
                    strand="+",
                    transcript_id="T1",
                    gene_id="G1",
                    gene_biotype="protein_coding",
                    transcript_support_level="1",
                ),
            ]
        )
        bam = self.make_bam(
            chroms=[("1", 1000)],
            segments=[
                dict(
                    qname="r1",
                    pos=150,
                    pnext=350,
                    cigar=[(0, 75)],
                    is_read1=True,
                    is_reverse=True,
                    **self.pkw,
                ),
                dict(
                    qname="r1",
                    pos=350,
                    pnext=150,
                    cigar=[(0, 75)],
                    is_read2=True,
                    mate_is_reverse=True,
                    **self.pkw,
                ),
            ],
        )

        ann1 = self.get_filename(extension="gtf")
        ann2 = self.get_filename(extension="gtf")
        result = self.get_filename(extension="gtf")

        make_conserved_gtf.process_gtf(gtf, ann1, source="ENSEMBL")
        make_conserved_intron_gtf.process_gtf(ann1, ann2)
        intronCounter_v2_stranded.intronCounter(bam, ann2, result)

        self.assertEqual(
            self.tsv_to_list(result, columns=[3, 4, 9, 10, 11]),
            [["200", "299", "1", "1", "5000000.0"]],
        )
