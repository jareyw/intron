##########
Change Log
##########

All notable changes to this project are documented in this file.


==========
Unreleased
==========

Changed
-------

Added
-----
- Add tox testing infrastructure. Tox enables us to run all tests with one
  single command: ``tox``. Tox runs:

    - unit tests
    - linter checks (black, pylint, pycodestyle, pydocstyle, isort)
    - checks that docs can be sucessfully built
    - checks that code is correctly packaged

- Add utilty classes and methods that help with generation of syntetic files
  used in tests

Fixed
-----
