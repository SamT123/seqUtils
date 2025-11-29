
<!-- README.md is generated from README.Rmd. Please edit that file -->

# seqUtils

<!-- badges: start -->

[![R-CMD-check](https://github.com/SamT123/seqUtils/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SamT123/seqUtils/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

Tools for working with biological sequences in R.

## Overview

`seqUtils` provides tools for common sequence analysis tasks (read &
write FASTA files, translate/align/compare sequences, etc.).

## Installation

You can install the development version of seqUtils from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("SamT123/seqUtils")
```

## Requirements

- MAFFT (for alignment)
- cmaple (for tree building)

## Example

``` r
library(seqUtils)

# Read FASTA file
sequences <- fast_fasta("sequences.fasta")

# Translate to amino acids
aa_sequences <- translate(sequences)

# Find substitutions compared to a reference
substitutions <- get_substitutions(
  reference_seq,
  query_seqs,
  exclude = c("X", "N")
)
```
