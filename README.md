# Binning and Evaluation Utilities

## Overview

This module contains the engines for the [MetaGenome Binning](https://github.com/BV-BRC/p3_binning) service and the
genome quality processor of the [Genome Annotation](https://github.com/BV-BRC/p3_genome_annotation) service.

## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

1. The [bins_generate.pl script](https://github.com/BV-BRC/p3_code/blob/master/scripts/bins_generate.pl) performs the
binning for bacteria and archaea.
2. The [vbins_generate.pl script](https://github.com/BV-BRC/p3_code/blob/master/scripts/vbins_generate.pl) performs the
binning for viruses.
3. The [BinningReports.pm module](https://github.com/BV-BRC/p3_code/blob/master/lib/BinningReports.pm) is used to
produce the binning reports.
4. The [p3x-eval-gto.pl script](https://github.com/BV-BRC/p3_code/blob/master/scripts/p3x-eval-gto.pl) is used to evaluate
the quality of the bacterial and archael bins.
5. The [p3x-improve-gto.pl script](https://github.com/BV-BRC/p3_code/blob/master/scripts/p3x-improve-gto.pl) removes contigs
from bins that are believed to be sources of contamination.
6. The [p3x-process-checkv.pl script](https://github.com/BV-BRC/p3_code/blob/master/scripts/p3x-process-checkv.pl) looks at
the results of the virus bin analysis to compute the actual bins.
7. The [eval_matrix.py script](https://github.com/BV-BRC/p3_code/blob/master/scripts/eval_matrix.py) uses RandomForest classifiers to
perform consistency checks on a genome's annotations.
