# SequenceBouncer 2.0

**SequenceBouncer** is a high-performance Python tool designed to identify and remove outlier or aberrant entries from multiple sequence alignments (MSAs). Outliers in MSAs can severely compromise downstream phylogenetic reconstruction and evolutionary analyses. SequenceBouncer uses column-based Shannon entropy combined with Hamming distance to fast, reproducibly, and context-responsively clean alignments.

Version 2.0 introduces **GPU (PyTorch CUDA) acceleration**, **multi-core CPU parallelism**, **memory-aware device fallbacks**, and **real-time progress bars**, achieving up to **1000x+ speedups** over the original implementation.

---

## Key Features in Version 2.0

*   **1000x+ Speedups:** Replaced nested loops with vectorized linear algebra (one-hot matrix multiplication), reducing calculations that originally took days to minutes or seconds.
*   **GPU Acceleration:** Offloads heavy matrix multiplications to NVIDIA Tensor Cores using PyTorch CUDA.
*   **Multi-Threaded CPU Parallelism:** Automatically utilizes optimized CPU BLAS backends (MKL, OpenBLAS) for parallelized CPU calculations.
*   **Memory-Aware Fallback:** Checks available GPU VRAM and CPU RAM before execution. If a full analysis exceeds GPU memory, it automatically falls back to CPU RAM rather than dropping to sampling mode, and warns the user.
*   **Thread & Device Control:** Allows users to explicitly limit CPU cores (`-c`/`--cores`) and select hardware targets (`-d`/`--device`).
*   **Real-time Progress Bars:** Displays visual status indicators using `tqdm` for sequence loading and sampling trials.
*   **Clean Output Directory Structure:** Writes all FASTA, CSV tables, logs, and PDF figures cleanly to an `output/` directory.

---

## Requirements & Installation

SequenceBouncer 2.0 is designed for Python 3.12+ and runs on Linux, macOS, and Windows.

We provide a pre-configured Conda environment file `environment.yml` to set up all dependencies (including PyTorch with CUDA 12.4 support).

### Setup environment:
```bash
# Clone the repository and navigate to the directory
cd SequenceBouncer

# Create the conda environment
conda env create -f environment.yml

# Activate the environment
conda activate sequencebouncer
```

---

## Usage

```bash
python SequenceBouncer.py -i <input_fasta> [options]
```

### Required Arguments:
*   `-i INPUT_FILE`, `--input_file INPUT_FILE`: Input multiple sequence alignment file in FASTA format.

### Optional Arguments:
*   `-o OUTPUT_FILE`, `--output_file OUTPUT_FILE`: Output filestem [do not include extensions] (default is `<input_filestem>`).
*   `-g GAP_PERCENT_CUT`, `--gap_percent_cut GAP_PERCENT_CUT`: Columns with a greater fraction of gaps than this percentage will be ignored in calculations (default is `2.0`).
*   `-k IQR_COEFFICIENT`, `--IQR_coefficient IQR_COEFFICIENT`: Coefficient multiplied by the interquartile range (IQR) to define an outlier cutoff (default is `1.0`).
*   `-n SUBSAMPLE_SIZE`, `--subsample_size SUBSAMPLE_SIZE`: Subsample size for sampling-based analysis (default is `0` for full alignment analysis). If GPU/CPU memory is insufficient, the system warns the user and automatically drops to CPU or prompts for sampling.
*   `-t TRIALS`, `--trials TRIALS`: Number of trials each sequence is sampled (default is `1` for full analysis, or `10` for sampling-based analysis).
*   `-s STRINGENCY`, `--stringency STRINGENCY`: Stringency for sampling-based filtering:
    *   `1`: Minimal stringency (retain sequences unless positive in all trials).
    *   `2`: Moderate stringency (retain sequences positive in $\le 50\%$ of trials).
    *   `3`: Maximum stringency (retain sequences positive in $0\%$ of trials).
*   `-r RANDOM_SEED`, `--random_seed RANDOM_SEED`: Integer seed for reproducible sampling-based shuffles.
*   `-c CORES`, `--cores CORES`: Maximum CPU threads/cores to use for CPU calculations. Defaults to **Available Cores - 4** to keep the host system responsive.
*   `-d DEVICE`, `--device DEVICE`: Choose calculation device:
    *   `auto` (default): Runs on GPU if available and memory fits, falling back to CPU if GPU memory is insufficient.
    *   `cpu`: Forces CPU calculation (using multi-threaded PyTorch CPU or NumPy).
    *   `cuda` or `gpu`: Forces GPU calculation, exiting with an error if CUDA is unavailable or VRAM is insufficient.
*   `-u`, `--unweighted`: Disable position-based inverse-entropy weighting and use standard unweighted Hamming distance (matching the original paper's behavior).
*   `-z`, `--gzip`: Force output clean and rejected FASTA files to be written in compressed gzip format (`.fasta.gz`), even if the input was uncompressed.
*   `-e ENTROPY_SMOOTHING`, `--entropy_smoothing ENTROPY_SMOOTHING`: Smoothing factor (pseudocount) added to Shannon entropy calculation to prevent division-by-zero and control position weighting strength (default is `0.01`). Higher values (e.g. `0.1` or `0.2`) soften the penalty assessed on conserved sites, protecting unusual but functional sequences from false rejection.

---

## Example Commands

### 1. Moderate Size Alignment (Full Analysis)
Runs a full pairwise alignment comparison using the CPU or GPU (device selected automatically):
```bash
python SequenceBouncer.py -i 1b_mammalian_mtDNA_input_alignment.fasta -k 1.75 -g 2.0
```

### 2. Large Alignment (Forced CPU Full Analysis)
Forces full analysis using exactly 8 CPU cores:
```bash
python SequenceBouncer.py -i zenodo_data/2b_bird_COI-5P_barcodes_input_alignment.fasta -d cpu -c 8
```

### 3. Extremely Large Alignment (GPU Sampling Analysis)
Runs a sampling-based analysis on the GPU with 10 trials and a subsample size of 100:
```bash
python SequenceBouncer.py -i zenodo_data/3b_COVID-19_input_alignment.fasta -n 100 -t 10 -s 2
```

---

## Output Files

All output files are saved in the `output/` directory:
*   `output/<filestem>_output_clean.fasta`: Contains the clean alignment with outlier sequences removed.
*   `output/<filestem>_output_rejected.fasta`: Contains the rejected outlier sequences.
*   `output/<filestem>_output_analysis.csv`: Tabular details showing total trials, positive outlier calls, and retention status.
*   `output/<filestem>_output.log`: Command logs and execution runtimes.
*   `output/<filestem>_gap_plot.pdf`: Plot showing gap fraction distribution across columns.
*   `output/<filestem>_median_plot.pdf` (Full analysis only): Plot showing sequence pairwise comparison medians and outlier cutoff.
*   `output/<filestem>_sampling_trials_plot.pdf` (Sampling mode only): Distribution of positive outlier calls across trials.

---

## Citation & Funding

If you use SequenceBouncer in your research, please cite:

> C.D. Dunn. (2020). SequenceBouncer: A method to remove outlier entries from a multiple sequence alignment. *bioRxiv* pre-print at: https://www.biorxiv.org/content/10.1101/2020.11.24.395459v2

### Funding & Acknowledgements
Development of this software was supported by funding from the **Sigrid Jusélius Foundation**. The author thanks Ed Davis (Oregon State University) for helpful suggestions.
