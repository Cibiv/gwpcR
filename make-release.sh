#!/bin/bash

set -e
set -o pipefail

ver=$1
if [ "$ver" == "" ]; then
	echo "Usage: $0 version" >&2
	exit 1
fi

# Update NAMESPACE and man/*.Rd files. If this changes anything, the check for uncomitted
# changes below will fire.
Rscript -e 'roxygen2::roxygenize()'

if [ "$(git symbolic-ref --short HEAD)" != "master" ]; then
	echo "Currently checkout out branch must be 'master'" >&2
	exit 1
fi

if [ $(git status --porcelain | wc -l) != 0 ]; then
	echo "Working copy contaisn uncommitted changes" >&2
	exit 1
fi

if [ $(git ls-remote --tags origin v$ver | wc -l) != 0 ]; then
	echo "Version $ver already released (remote tag v$ver already exists)" >&2
	exit 1
fi

# Run tests
Rscript -e 'devtools::test()'

echo "Updating DESCRIPTION" >&2
sed -i.bak 's/^Version: \(.*\)$/Version: '"$ver"'/' DESCRIPTION
rm DESCRIPTION.bak

echo "Comitting change" >&2
git add DESCRIPTION
git commit -m "Incremented version to $ver"

echo "Tagging as v$ver and latest" >&2
git tag -f v$ver

echo "Pushing to origin" >&2
git push origin master v$ver

echo "Updating 'latest' on origin" >&2
git push -f origin v$ver:refs/tags/latest-release
