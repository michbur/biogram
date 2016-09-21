[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/biogram)](https://cran.r-project.org/package=biogram)
[![Downloads](http://cranlogs.r-pkg.org/badges/biogram)](https://cran.r-project.org/package=biogram)
[![Build Status](https://api.travis-ci.org/michbur/biogram.png)](https://travis-ci.org/michbur/biogram)
[![codecov.io](https://codecov.io/github/michbur/biogram/coverage.svg?branch=master)](https://codecov.io/github/michbur/biogram?branch=master) 


biogram package
------------

This package contains tools for extraction and analysis of various
n-grams (sequences of n items) derived from biological sequences (proteins
or nucleic acids). To deal with the curse of dimensionality of the n-grams,
biogram uses Quick Permutation Test (QuiPT) for fast feature filtering.

Installation
------------

biogram is available [on CRAN](https://cran.r-project.org/package=biogram), so installation is as simple as:

```
install.packages("biogram")
```

You can install the latest development version of the code using the `devtools` R package.

```
# Install devtools, if you haven't already.
install.packages("devtools")

library(devtools)
install_github("michbur/biogram")
```

For citation type:

```
citation("biogram")
```

or use:
Michal Burdukiewicz, Piotr Sobczyk and Chris Lauber (2016). biogram: N-Gram Analysis of Biological Sequences. R package version 1.3. https://cran.r-project.org/package=biogram
