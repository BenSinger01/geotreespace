# *geotreespace*: an R Package for the Statistical Measure of Phylogeographic Incopatibilities

*geotreespace* implements new mathematical measures of phylogeographic concordance or incompatibilities,:
the Pairwise Most-Recent Common Ancestor (MRCA) incompatibility $I_{mrca}$ and the Maximum Agreement Sub-Phylogeography (MASPG) $I_{maspg}$.
These emtrics are aimed at detecting differences between geographical histories in terms of
distances between phylogeographies. These incompatibility measures can be applied
to any phylogeography, and more generally to any phylogeny where each tip has been assigned
a geographical discrete location or any other continuous/discrete 'trait' independent of the sequence.
A non-zero value of an incompatibility measure between two phylogeographies implies that they are incompatible,
and could be evidence of recombination.


## Installing *geotreespace*

To install the development version from github:

    library(devtools)
    install_github("BenSinger01/geotreespace")

Then, to load the package, use:
    library("geotreespace")


## Content overview

The main functions implemented in *geotreespace* are:

-  **`mrcaDist`**: comparison of two phylogeographic trees using the MRCA incompatibility measure

-  **`multiMrcaDist`**: comparison of a list of phylogeographic trees using the MRCA incompatibility measure

-  **`maspgDist`**: comparison of two phylogeographic trees using the MASPG incompatibility metric

-  **`multiMaspgDist`**: Comparison of a list of phylogeographic trees using the MASPG incompatibility metric

Other functions are central to the computations of distances between
trees:

-  **`mrcaVec`**: function to find the location of the most recent common ancestor (MRCA) of each pair of tips in a tree


## Reporting bugs

Bugs can be reported using the package’s
[issue system](https://github.com/BenSinger01/geotreespace/issues).


## Authors

Authors:

- [Benjamin John Singer](https://twitter.com/_bensinger?lang=en-GB)

- [Antonello Di Nardo](https://www.pirbright.ac.uk/users/dr-antonello-di-nardo)

- [Luca Ferretti](https://sites.google.com/view/lucaferretti)
