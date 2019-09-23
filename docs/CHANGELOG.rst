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
- Rewrite intronCounter script into module and expose it as entry point.
  This enables dual usage::

    - It can be imported and used as any ordinary Python function
    - It can be run from command line with ``intron-counter`` command

Fixed
-----
