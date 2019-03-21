.. _development-guidelines:

Development guidelines
======================

If you want to contribute code to Ocellaris then, first of all, Thank You! 🙌
Please break your changes into reasonably sized pull requests and follow the
rather standard and simple guidelines below. If you are afraid of doing
something wrong, don't let that stop you! As long as you have made an effort to
follow the expected steps then you will get feedback and we will get any
problems fixed. There is no guarantee that your pull request will be merged,
but not submitting 1000 changed files with ugly code and no tests will help
your case significantly 😉

.. contents:: Steps to get your code merged
    :local:

Please keep in mind that the person/persons that can review your pull request
may not get to it right away. In the worst case scenario they are on a month
long tropical Holiday away from email (best case scenario! ☀️🌴🍹). Do not get
discouraged if we do not get to reviewing your changes right away!


Git workflow
------------

Ocellaris is developed in git using the standard `feature branch workflow
<https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow>`_.
Only documentation commits should go directly to the master branch, all other
changes should be made in their own branch (branched off master) and then
merged to master when they are ready through a pull request.

The following are examples of good branch names:

* tormod/implement-foobar-bcs
* fred/pressure-experiments
* fredrickson/issue-123-fix_barfoo
* amy/issue-246

Prefixing the branch name with your own name helps others know who is the main
responsible for each branch, which can be useful if more people want to
contribute to the Ocellaris development.

You will get write access to the Ocellaris repositories if you ask nicely and
have contributed a couple of pull requests. You will not be given access to
merge to the master branch at once, though.


Code style
----------

Ocellaris' code follows the standard Python PEP8 code style, more or less.
Instead of describing the code style in detail it should be sufficient to
run

    `flake8 <https://pypi.org/project/flake8/>`_ ``--max-line-length=100
    --ignore=E203``  *YOUR_CODE.py*

which will report any errors. Normally you do not run this yourself, but
configure your editor to run flake8 with these swiches for you.

Since you probably have better things to do than wory about the formatting of
your code you can run

    `black <https://pypi.org/project/black/>`_ ``--line-length=100
    --skip-string-normalization``  *YOUR_CODE.py*

Black is a code formatter and can be configured to run inside your editor every
time you save. That way your code can look bad when you write it, but, as long
as it is valid Python 3 code, black will fix most issues after you press
``Ctrl+S``. Black will not fix or format syntaxtically incorrect code, and it
will not fix all errors raised by flake8, just most of them.


Automated tests
---------------

Testing numerical code can be challenging. When a large amount of global state
is involved (all degrees of freedom in the problem, ``simulation.data[*]``,
...) it can easily happen that a change in one place of the code can have
unintended consequences somewhere else. Please try to include unit tests if
your code has functions or classes that do something testable without running a
full simulation (should be most code) and consider using MMS tests or short
demos with expected results to show that the changes perform as expected in a
full simulation as well.

You must make sure that existing tests run before merging your code. The
automatic testing will tell you if you broke something, but it is easier to
debug the tests locally. See the CircleCI configuration in the code repository
for how to run the tests.

Please resist the temptation to include new test that take forever to complete.
For the sake of everyone in the future having to wait for that test, please try
to write tests that run relatively fast. Unit tests should mostly be
instantaneous and the integration tests that run a full simulation should not
take more than 30 seconds. Less is better! Lets all try to keep the tests
efficient wile still ensuring that the code is correct. Need a mesh? Maybe a
small 2D mesh is sufficient for the test?

Ocellaris uses `pytest <https://pytest.org/>`_ to run both unit and integration
tests. If you are not familiar with unit tests then just google what they are
and you will find about a billion articles extolling their virtues. Think of
the test you write as helpers to protect you against future programmers
breaking your code because they do not understand properly how it works. Often
that future person is yourself two weeks after ... 🤣

.. _label-running-tests:

**Running the tests**

The exact procedure for running the tests can be found in `config.yml
<https://bitbucket.org/ocellarisproject/ocellaris/src/master/.circleci/config.yml>`_.
Currently this is the following commands, run them in the repository base directory:

.. code-block:: bash

    # Install test dependencies
    python3 -m pip install pytest pytest-timeout pytest-instafail

    # Run unit tests
    python3 -m pytest -v tests/ --instafail --timeout=300 --durations=10

    # Run unit tests with MPI
    mpirun -n 3 python3 -m pytest -v tests/ --maxfail=1 --timeout=300 --durations=10
    
    # Run regression tests
    python3 -m pytest -v cases/regression_tests.py --instafail --timeout=300 --durations=0
    
    # Run demo tests
    python3 -m pytest -v demos/ --instafail --timeout=300 --durations=0
