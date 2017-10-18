#!/bin/bash
# Build documentation and upload to bitbucket pages
# Runs as a part of the automated testing on CircleCI


# Script setup
set -euxo pipefail
REPO=trlandet@bitbucket.org:trlandet/trlandet.bitbucket.io.git
CMESSAGE="CircleCI Ocellaris doc update"
cd documentation/


# Install and run Sphinx
pip3 install sphinx --force --user
export PATH=~/.local/bin:$PATH
make html


# Configure git
git config --global user.email "nobody@example.com"
git config --global user.name "CircleCI Ocellaris Doc Builder"


# Clone and update the webpage git repo
git clone --depth 1 $REPO webpage
rm -r webpage/ocellaris
mv _build/html webpage/ocellaris


# Add doc changes to repo
cd webpage/
git add -A


# Exit if there are no changes
if git diff-index --quiet HEAD --; then
  exit 0
fi


# Amend previous commit if it was a CI commit, else
# make a new commit (to avoid filling up the repo
# with autogenerated commits
if [[ "$(git log -1 --pretty=%B)" == *"$CMESSAGE"* ]]; then
  git commit --amend -m "$CMESSAGE"
  git push -f
else
  git commit -m "$CMESSAGE"
  git push
fi