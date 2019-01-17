.. _development-guidelines:

Development guidelines
======================

Ocellaris is developed in git using the standard `feature branch workflow
<https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow>`_.
Only documentation commits should go directly to the master branch, all other
changes should be made in their own branch (branched off master) and then
merged to master when they are ready through a pull request.

The following are examples of good branch names:

* tormod/implement-foobar-bcs
* tormod/pressure-experiments
* tormod/issue-123-fix_barfoo
* tormod/issue-123

Prefixing the branch name with your own name helps others know who is the main
responsible for each branch, which can be useful if more people want to
contribute to the Ocellaris development.
