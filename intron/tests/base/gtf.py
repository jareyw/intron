"""Test utilities for GTF files."""
import pybedtools as pbt


class GtfTestCaseMixin:
    """Mixin for manipulation of GTF files in tests."""

    def create_gtf_attrs(self, **kwargs):
        """Create 9th column in GTF file - attributes."""
        attrs = ""
        for name, value in kwargs.items():
            attrs += '{} "{}"; '.format(name, value)

        attrs = attrs.strip()
        if attrs:
            return attrs
        else:
            return "."

    def make_gtf_interval(
        self,
        seqname="chr1",
        source="UCSC",
        feature="exon",
        start=1,
        end=1000,
        score=1,
        strand="+",
        frame=".",
        **kwargs,
    ):
        """Create GTF interval."""
        attrs = self.create_gtf_attrs(**kwargs)
        # GTF is 1 based so start+gtf=start_python + 1
        # GTF does include end coordinate, but Python does not, so end_gtf == end_python
        data = [seqname, source, feature, start + 1, end, score, strand, frame, attrs]
        return pbt.create_interval_from_list(data)

    def make_gtf(self, intervals):
        """Create GTF file from a list of lists.

        Data should be a list of pbt.Interval objects.
        """
        fname = self.get_filename(extension="gtf")
        pbt.BedTool(interval for interval in intervals).saveas(fname)
        return fname
