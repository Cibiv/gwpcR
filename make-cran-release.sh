#!/bin/bash

set -e
set -o pipefail

test -e ".cran-pkg" && (echo "directory .cran-pkg already exists"; exit 1)

mkdir .cran-pkg

cp -r 	DESCRIPTION \
	NAMESPACE \
	COPYING \
	R \
	src \
	man \
	inst \
	tests \
	.cran-pkg

R CMD build .cran-pkg

rm -r .cran-pkg

echo "DONE"


