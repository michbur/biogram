biogram package
------------

This package contains tools for extraction and analysis of various
n-grams (sequences of n items) derived from biological sequences (proteins
or nucleic acids). To deal with the curse of dimensionality of the n-grams,
biogram uses Quick Permutation Test (QuiPT) for fast feature filtering.

Installation
------------

biogram is available [on CRAN](http://cran.r-project.org/web/packages/biogram/), so installation is as simple as:

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
Michal Burdukiewicz, Piotr Sobczyk and Chris Lauber (2015). biogram: N-Gram Analysis of Biological Sequences. R package version 1.2. http://CRAN.R-project.org/package=biogram
