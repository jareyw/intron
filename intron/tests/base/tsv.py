"""Test utilities for .tsv files."""

import gzip


class TsvTestCaseMixin:
    """Mixin for manipulation of .tsv files in tests."""

    def tsv_to_list(self, fname, fields_separator="\t", columns=None):
        """Read .tsv file to a nested list."""
        if columns and not isinstance(columns, list):
            columns = [columns]

        open_function = open
        if fname.endswith(".gz"):
            open_function = gzip.open

        content = []
        with open_function(fname, "rt") as handle:
            for line in handle:
                line = line.strip().split(fields_separator)

                if columns:
                    line = [item for i, item in enumerate(line) if i in columns]

                content.append(line)
        return content
