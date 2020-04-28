##########
Change Log
##########

All notable changes to this project are documented in this file.


==================
0.1.0 - 2020-04-28
==================

Changed
-------
- Change directory structure
- Change testing configuration

Added
-----
- Support analysis with ENSEMBL annotation

Fixed
-----
- Remove too long genes - the UCSC ones that could be different genes
  with same id


==================
0.0.1 - 2020-01-31
==================

Changed
-------
- Port the following scripts to Python3, put them in the main intron folder
  and make some minor modifications::

  - make_conserved_gtf.py
  - make_conserved_intron_gtf.py
  - intronCounter_v2_stranded.py
- Remove legacy folder

Added
-----
- Add tox testing infrastructure. Tox enables us to run all tests with one
  single command: ``tox``. Tox runs:

    - unit tests
    - linter checks (black, pylint, pycodestyle, pydocstyle, isort)
    - checks that docs can be sucessfully built
    - checks that code is correctly packaged
- Add semaphore configuration
- Add some basic tests
