SequenceBouncer
Version: 1.3

----

Author:

Cory Dunn
Contact: cory.david.dunn@gmail.com

----

License:

GPLv3

----

Phylogenetic analyses can take advantage of multiple sequence alignments as input. These alignments typically consist of homologous nucleic acid or protein sequences, and the inclusion of outlier or aberrant sequences can compromise downstream analyses. SequenceBouncer uses the Shannon entropy values of alignment columns to identify, then quickly and reproducibly remove, outlier entries in a manner responsive to overall alignment context.

----

Please cite: 

C.D. Dunn. (2020). SequenceBouncer: A method to remove outlier entries from a multiple sequence alignment. bioRxiv pre-print at: https://www.biorxiv.org/content/10.1101/2020.11.24.395459v2

----

Acknowledgements:

Funding received from the Sigrid Jusélius Foundation contributed to the development of this software. The author thanks Ed Davis, Oregon State University, for helpful suggestions.

----

Requirements:

SequenceBouncer is implemented in Python 3 (current version tested under Python version 3.12.8) 

Dependencies: 

Biopython (tested under version 1.85),
Pandas (tested under version 2.2.3),
Numpy (tested under version 1.26.4),
Matplotlib (tested under version 3.10.0)

----

Usage:

SequenceBouncer.py [-h] -i INPUT_FILE [-o OUTPUT_FILE] [-g GAP_PERCENT_CUT]
                                             [-k IQR_COEFFICIENT] [-n SUBSAMPLE_SIZE] [-t TRIALS]
                                             [-s STRINGENCY] [-r RANDOM_SEED]

Required arguments:
  -i INPUT_FILE, --input_file INPUT_FILE
                        Input multiple sequence alignment file in FASTA format.

Optional arguments:

  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output filename [do not include extensions] (default will be 'input_file.ext').

  -g GAP_PERCENT_CUT, --gap_percent_cut GAP_PERCENT_CUT
                        For columns with a greater fraction of gaps than the selected value, expressed in
                        percent, data will be ignored in calculations (default is 2).

  -k IQR_COEFFICIENT, --IQR_coefficient IQR_COEFFICIENT
                        Coefficient multiplied by the interquartile range that helps to define an outlier
                        sequence (default is 1.0).

  -n SUBSAMPLE_SIZE, --subsample_size SUBSAMPLE_SIZE
                        |> Available for large alignments | The size of a single sample taken from the
                        full dataset (default is entire alignment, but try a subsample size of 50 or 100
                        for large alignments).

  -t TRIALS, --trials TRIALS
                        |> Available for large alignments | Number of times each sequence is sampled and
                        tested (default is to examine all sequences in one single trial, but 5 or 10
                        trials may work well when subsamples are taken from large alignments).

  -s STRINGENCY, --stringency STRINGENCY
                        |> Available for large alignments | 1: Minimal stringency 2: Moderate stringency
                        3: Maximum stringency (default is moderate stringency).

  -r RANDOM_SEED, --random_seed RANDOM_SEED
                        Random seed (integer) to be used during a sampling-based approach (default is
                        that the seed is randomly selected). The user can use this seed to obtain
                        reproducible output and should note it in their publications.

---

This could be the starting point for further exploration of parameters:

i) If the alignment is of moderate size: 

python SequenceBouncer.py -i <input alignment>

ii) If the alignment is of substantial size:

python SequenceBouncer.py -i <input alignment> -n 100 -t 10 -s 2

---

Example files are found at: 

http://doi.org/10.5281/zenodo.4285789
                        
