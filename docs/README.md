# Guenomu documentation

Guenomu relates gene-family histories to a shared species history using Bayesian tree reconciliation. Input consists of gene-tree collections and species labels; output includes sampled species trees, gene trees and reconciliation parameters.

Start with [installation](installation.md), then copy the supplied example data into a working directory and follow the [output guide](output.md). The [model paper](papers/paper001_guenomu.md) and [software chapter](papers/paper002_guenomu.md) describe the published statistical method and practical workflow. They are preserved as historical publications; the installation and output guides describe current behavior.

The C library includes tree handling, reconciliation, random-number generation and numerical utilities. Its interface remains experimental; see [development](development.md).
