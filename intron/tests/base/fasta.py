"""Test utilities for FASTA files."""
import os

import numpy

BASES = ["A", "C", "G", "T"]


class FastaTestCaseMixin:
    """Mixin for manipulation of FASTA files in tests."""

    @staticmethod
    def make_fasta_sequence(size, include_n=False, rnd_seed=None):
        """Make random DNA segment of length `size`."""
        numpy.random.seed(rnd_seed)
        bases = BASES + ["N"] if include_n else BASES
        return "".join(numpy.random.choice(bases, size))

    def make_fasta(self, sequences=None, headers=None, num_sequences=10, seq_len=80, rnd_seed=None):
        """Make FASTA file."""
        if sequences and headers is None:
            num_sequences = len(sequences)
        if headers and sequences is None:
            num_sequences = len(headers)

        numpy.random.seed(rnd_seed)
        if sequences is None:
            random_seeds = numpy.random.randint(10 ** 5, size=num_sequences)
            sequences = [self.make_fasta_sequence(seq_len, rnd_seed=rnd) for rnd in random_seeds]
        if headers is None:
            headers = ["{}".format(i + 1) for i in range(num_sequences)]

        out_file = self.get_filename(extension="fasta")
        with open(out_file, "wt") as ofile:
            for header, seq in zip(headers, sequences):
                ofile.write(">" + header + "\n")
                ofile.write(seq + "\n")

        return os.path.abspath(out_file)
