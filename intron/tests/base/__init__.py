"""Testing utilities."""
import os
import tempfile
import unittest

from .bam import BamTestCaseMixin
from .fasta import FastaTestCaseMixin
from .gtf import GtfTestCaseMixin
from .tsv import TsvTestCaseMixin


class IntronTestCase(BamTestCaseMixin, FastaTestCaseMixin, GtfTestCaseMixin, TsvTestCaseMixin, unittest.TestCase):
    """Base class for intron TestCases."""

    @staticmethod
    def get_filename(extension=None, directory=None):
        """Get availbale filename."""
        name = next(tempfile._get_candidate_names())

        if not directory:
            directory = tempfile._get_default_tempdir()
        if extension:
            name += "." + extension

        return os.path.join(directory, name)

    @staticmethod
    def get_tmp_dir():
        """Return a temporary directory."""
        return tempfile.mkdtemp()
