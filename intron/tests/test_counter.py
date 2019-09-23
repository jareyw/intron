"""Intron counter tests."""

from intron.counter import count_introns

from .base import IntronTestCase


class CounterTest(IntronTestCase):
    """Base class for intron TestCases."""

    def setUp(self):
        """Set up testing files."""
        self.outfile = self.get_filename("tsv.gz")

        self.gtf = self.make_gtf(
            [
                # Positive strand
                self.make_gtf_interval(feature="intron", start=100, end=200, strand="+"),
                # Negative strand
                self.make_gtf_interval(feature="intron", start=400, end=500, strand="-"),
                # Short interval
                self.make_gtf_interval(feature="intron", start=800, end=810, strand="+"),
            ]
        )

        # Paired kwargs:
        self.pkw = {"is_paired": True, "is_proper_pair": True}

    def test_basic(self):
        """Basic test."""
        bam = self.make_bam(
            chroms=[("chr1", 1000), ("chr2", 1000)],
            segments=[
                dict(qname="r001", pos=150, pnext=350, cigar=[(0, 100)], is_read1=True, **self.pkw),
                dict(qname="r001", pos=350, pnext=150, cigar=[(0, 100)], is_read2=True, **self.pkw),
            ],
        )

        count_introns(bam, self.gtf, self.outfile)
        self.assertEqual(
            self.tsv_to_list(self.outfile, columns=[0, 3, 4, 9, 10, 11]),
            [
                ["chr1", "101", "200", "1", "0", "5000000.0"],
                ["chr1", "401", "500", "0", "0", "0.0"],
                ["chr1", "801", "810", "0", "0", "0.0"],
            ],
        )

    def test_short_intron(self):
        """Test that too short introns have zero score."""
        bam = self.make_bam(
            chroms=[("chr1", 1000), ("chr2", 1000)],
            segments=[
                dict(qname="r001", pos=750, pnext=850, cigar=[(0, 100)], is_read1=True, **self.pkw),
                dict(qname="r001", pos=850, pnext=750, cigar=[(0, 100)], is_read2=True, **self.pkw),
            ],
        )

        count_introns(bam, self.gtf, self.outfile)
        self.assertEqual(
            self.tsv_to_list(self.outfile, columns=[0, 3, 4, 9, 10, 11]),
            [
                ["chr1", "101", "200", "0", "0", "0.0"],
                ["chr1", "401", "500", "0", "0", "0.0"],
                ["chr1", "801", "810", "0", "0", "0.0"],
            ],
        )

    def test_span(self):
        """Test difference between spanning read and read aligning to intron."""
        bam = self.make_bam(
            chroms=[("chr1", 1000), ("chr2", 1000)],
            segments=[
                dict(qname="r001", pos=150, pnext=350, cigar=[(3, 50), (0, 50)], is_read1=True, **self.pkw),
                dict(qname="r001", pos=350, pnext=150, cigar=[(0, 100)], is_read2=True, **self.pkw),
            ],
        )

        count_introns(bam, self.gtf, self.outfile)
        self.assertEqual(
            self.tsv_to_list(self.outfile, columns=[0, 3, 4, 9, 10, 11]),
            [
                ["chr1", "101", "200", "0", "1", "0.0"],
                ["chr1", "401", "500", "0", "0", "0.0"],
                ["chr1", "801", "810", "0", "0", "0.0"],
            ],
        )
