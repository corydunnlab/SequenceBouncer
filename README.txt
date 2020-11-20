SequenceBouncer
----
Author:
Cory Dunn
Institute of Biotechnology
University of Helsinki
Email: cory.dunn@helsinki.fi
----
License:
GPLv3
----
Acknowledgements:
Funding received from the Sigrid Jusélius Foundation contributed to the development of this software.
Code utilized for Shannon entropy calculation, written by Joe R. J. Healey, University of Warwick (https://gist.github.com/jrjhealey/) was modified and used under the GPLv3 license.
----
Requirements:
SequenceBouncer is implemented in Python 3 (tested under version 3.8.6) 
Dependencies: 
Biopython (tested under version 1.78),
Pandas (tested under version 1.1.3),
Numpy (tested under version 1.19.1)
----
Usage:
SequenceBouncer_v1_16.py [-h] -i INPUT_FILE [-o OUTPUT_FILE] [-s STRINGENCY] [-k IQR_COEFFICIENT] [-t TRIALS] [-n SUBSAMPLE_SIZE] [-g GAP_PERCENT_CUT]

Required arguments:
  -i INPUT_FILE, --input_file INPUT_FILE
                        Input file in FASTA format.

Optional arguments:
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output filename [do not include extensions] (Default will be 'input_file.ext').

  -s STRINGENCY, --stringency STRINGENCY
                        1: Minimal stringency. 2: Moderate stringency. 3: Maximum stringency. (Default is moderate stringency).


  -k IQR_COEFFICIENT, --IQR_coefficient IQR_COEFFICIENT
                        Coefficient multiplied by the interquartile range that helps to define an outlier sequence (Default is 1.0)..

  -t TRIALS, --trials TRIALS
                        Number of times each sequence is sampled and tested (Default is to examine all sequences in one single
                        trial).

  -n SUBSAMPLE_SIZE, --subsample_size SUBSAMPLE_SIZE
                        The size of a single sample taken from the full dataset (Default is entire alignment).

  -g GAP_PERCENT_CUT, --gap_percent_cut GAP_PERCENT_CUT
                        For columns with a greater fraction of gaps than the selected value (expressed in percent), data will be ignored in calculations (Default is 2).
                        