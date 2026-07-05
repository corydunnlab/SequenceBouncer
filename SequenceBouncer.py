#!/usr/bin/env python
# coding: utf-8 

# Funding received from the Sigrid Jusélius Foundation contributed to the development of this software.
# Author: Cory Dunn
# Author Email: cory.david.dunn@gmail.com
# Version: 2.0 (Memory Optimized)
version = '2.0-memopt'
# License: GPLv3

import sys
import os

# 1. Early parsing of cores argument to set environment variables before numpy/pytorch imports
cores = None
for i, arg in enumerate(sys.argv):
    if arg in ('-c', '--cores'):
        if i + 1 < len(sys.argv):
            try:
                cores = int(sys.argv[i+1])
            except ValueError:
                pass
        break

if cores is None:
    try:
        cores = max(1, os.cpu_count() - 4) if os.cpu_count() is not None else 1
    except Exception:
        cores = 1

# Set threads for BLAS libraries before numpy is loaded
os.environ["OMP_NUM_THREADS"] = str(cores)
os.environ["MKL_NUM_THREADS"] = str(cores)
os.environ["OPENBLAS_NUM_THREADS"] = str(cores)
os.environ["NUMEXPR_NUM_THREADS"] = str(cores)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(cores)

# Now safe to import pandas, numpy, torch, etc.
from Bio import SeqIO
import pandas as pd
import numpy as np
import time
import argparse
import gc
import random
import logging
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Try to import tqdm for progress bar
try:
    from tqdm import tqdm
    tqdm_available = True
except ImportError:
    tqdm_available = False

if not tqdm_available:
    def tqdm(iterable, *args, **kwargs):
        return iterable

# Try to import torch for GPU acceleration
try:
    import torch
    torch_available = True
    # Limit PyTorch CPU threads to the selected cores count
    torch.set_num_threads(cores)
except ImportError:
    torch_available = False

def get_torch_device():
    if torch_available:
        if torch.cuda.is_available():
            return torch.device("cuda")
        elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            return torch.device("mps")
    return None

def get_available_memory():
    gpu_mem = None
    device = get_torch_device()
    
    if device is not None:
        if device.type == "cuda":
            try:
                total_vram = torch.cuda.get_device_properties(device).total_memory
                gpu_mem = total_vram * 0.7
            except Exception:
                pass
        elif device.type == "mps":
            # Apple Silicon uses Unified Memory, so VRAM is shared with system RAM
            try:
                import psutil
                gpu_mem = psutil.virtual_memory().available * 0.7
            except ImportError:
                pass
            
    # Fallback to CPU RAM on Linux/macOS
    cpu_mem = 8 * 1024 * 1024 * 1024 # default fallback 8 GB
    if os.path.exists("/proc/meminfo"):
        try:
            with open("/proc/meminfo", "r") as f:
                for line in f:
                    if "MemAvailable" in line:
                        parts = line.split()
                        cpu_mem = int(parts[1]) * 1024 # convert KB to bytes
                        break
        except Exception:
            pass
    else:
        try:
            import psutil
            cpu_mem = psutil.virtual_memory().available
        except ImportError:
            pass
            
    return gpu_mem, cpu_mem

def compute_distances(table_sample_numpy, max_entropy_after_gaps, weights=None, force_cpu=False):
    N, L = table_sample_numpy.shape
    device = get_torch_device()
    
    # Try GPU via PyTorch if GPU (CUDA/MPS) is available, matrix is large, and not forced to CPU
    if device is not None and N >= 200 and not force_cpu:
        try:
            t = torch.from_numpy(table_sample_numpy).to(device)
            unique_vals = torch.unique(t)
            one_hot = (t.unsqueeze(-1) == unique_vals).float()
            if weights is not None:
                t_w = torch.from_numpy(weights).to(device).float()
                scaled_one_hot = one_hot * torch.sqrt(t_w).unsqueeze(-1)
                scaled_one_hot = scaled_one_hot.view(N, -1)
                matches = torch.matmul(scaled_one_hot, scaled_one_hot.T)
                total_weight = torch.sum(t_w)
                transformed_distances = total_weight - matches
            else:
                one_hot = one_hot.view(N, -1)
                matches = torch.matmul(one_hot, one_hot.T)
                mismatches = L - matches
                transformed_distances = max_entropy_after_gaps * mismatches
            return transformed_distances.cpu().numpy()
        except Exception:
            pass
            
    # Try PyTorch CPU if PyTorch is available but GPU is not (often faster than NumPy for large arrays)
    if torch_available and N >= 200:
        try:
            t = torch.from_numpy(table_sample_numpy)
            unique_vals = torch.unique(t)
            one_hot = (t.unsqueeze(-1) == unique_vals).float()
            if weights is not None:
                t_w = torch.from_numpy(weights).float()
                scaled_one_hot = one_hot * torch.sqrt(t_w).unsqueeze(-1)
                scaled_one_hot = scaled_one_hot.view(N, -1)
                matches = torch.matmul(scaled_one_hot, scaled_one_hot.T)
                total_weight = torch.sum(t_w)
                transformed_distances = total_weight - matches
            else:
                one_hot = one_hot.view(N, -1)
                matches = torch.matmul(one_hot, one_hot.T)
                mismatches = L - matches
                transformed_distances = max_entropy_after_gaps * mismatches
            return transformed_distances.numpy()
        except Exception:
            pass
            
    # Standard NumPy implementation (highly optimized via BLAS)
    unique_chars = np.unique(table_sample_numpy)
    if weights is not None:
        one_hot = (table_sample_numpy[:, :, None] == unique_chars[None, None, :])
        scaled_one_hot = one_hot * np.sqrt(weights)[None, :, None]
        scaled_one_hot = scaled_one_hot.reshape(N, -1).astype(np.float32)
        matches = np.dot(scaled_one_hot, scaled_one_hot.T)
        total_weight = np.sum(weights)
        transformed_distances = total_weight - matches
    else:
        one_hot = (table_sample_numpy[:, :, None] == unique_chars[None, None, :]).reshape(N, -1).astype(np.float32)
        matches = np.dot(one_hot, one_hot.T)
        mismatches = L - matches
        transformed_distances = max_entropy_after_gaps * mismatches
    return transformed_distances

if __name__ == "__main__" :
    default_cores = max(1, os.cpu_count() - 4) if os.cpu_count() is not None else 1
    
    ap = argparse.ArgumentParser()
    ap.add_argument('-i','--input_file',required=True,type=str,help='Input file in FASTA format.\n')
    ap.add_argument('-o','--output_file',required=False,type=str,default='X',help="Output filestem [do not include extensions] (default will be '<input_filestem>.ext').\n")
    ap.add_argument('-g','--gap_percent_cut',required=False,type=float,default=2.0,help='For columns with a greater fraction of gaps than the selected value, expressed in percent, data will be ignored in calculations (default is 2).\n')
    ap.add_argument('-k','--IQR_coefficient',required=False,type=float,default=1.0,help='Coefficient multiplied by the interquartile range that helps to define an outlier sequence (default is 1.0).\n')
    ap.add_argument('-n','--subsample_size',required=False,type=int,default=0,help='|> Available for large alignments | The size of a single sample taken from the full dataset (default is 0 for auto/full alignment, automatically falling back to CPU or sampling based on resources).\n')
    ap.add_argument('-t','--trials',required=False,type=int,default=0,help='|> Available for large alignments | Number of times each sequence is sampled and tested (default is 1 for full analysis, or dynamically set for sampling analysis).\n')
    ap.add_argument('-s','--stringency',required=False,type=int,default=2,help='|> Available for large alignments | 1: Minimal stringency 2: Moderate stringency 3: Maximum stringency (default is moderate stringency).\n')
    ap.add_argument('-r','--random_seed',required=False,type=int,default=random.randint(0,1000),help='Random seed (integer) to be used during a sampling-based approach (default is that the seed is randomly selected). The user can use this seed to obtain reproducible output and should note it in their publications. \n')
    ap.add_argument('-c','--cores',required=False,type=int,default=default_cores,help='Number of CPU threads to use for calculations (default is available number minus 4).\n')
    ap.add_argument('-d','--device',required=False,type=str,default='auto',choices=['auto','cpu','cuda','gpu','mps'],help='Device to use for calculations: auto, cpu, cuda, gpu, mps (default is auto).\n')
    ap.add_argument('-m','--mode',required=False,type=str,default='auto',choices=['auto','full','sample'],help='Analysis mode: auto (automatically choose based on resources), full (force full pairwise analysis), sample (force sampling-based analysis). Default is auto.\n')
    ap.add_argument('-u','--unweighted',action='store_true',help='Disable position-based weighting and use standard unweighted Hamming distance instead.\n')
    ap.add_argument('-z','--gzip',action='store_true',help='Force output clean and rejected FASTA files to be written in gzipped format (.fasta.gz).\n')
    ap.add_argument('-e','--entropy_smoothing',required=False,type=float,default=0.01,help='Smoothing factor (pseudocount) added to Shannon entropy calculation to prevent division-by-zero and control the strength of position-based weighting (default is 0.01). Higher values like 0.1 soften the penalty on conserved sites.\n')

    args = vars(ap.parse_args())
    input_sequence = args['input_file']
    stringency_flag = args['stringency']
    min_trials_for_each_sequence = args['trials']
    multiplier_on_interquartile_range = args['IQR_coefficient']
    number_in_small_test = args['subsample_size']
    gap_value_cutoff = args['gap_percent_cut']
    output_entry = args['output_file']
    seed = args['random_seed']
    cores = args['cores']
    device_flag = args['device'].lower()
    mode_flag = args['mode'].lower()
    weighted = not args['unweighted']
    gzip_flag = args['gzip']
    entropy_smoothing = args['entropy_smoothing']

    # Dynamic target resolution for 'gpu'
    target_device = get_torch_device()
    if device_flag == 'gpu':
        if target_device is not None:
            device_flag = target_device.type
        else:
            device_flag = 'cuda'  # default fallback

    # Set threads again to guarantee PyTorch limits inside main execution
    if torch_available:
        torch.set_num_threads(cores)

    # Setup directories
    WD = os.getcwd()
    output_directory = os.path.join(WD, 'output')
    os.makedirs(output_directory, exist_ok=True)

    if output_entry == 'X':
        # Find input basename without extensions like .fasta, .fasta.gz, .gz
        base = os.path.basename(input_sequence)
        if base.endswith('.gz'):
            base = base[:-3]
        for ext in ('.fasta', '.fa', '.fna', '.alignment', '.aln'):
            if base.lower().endswith(ext):
                base = base[:-len(ext)]
                break
        output_entry = base

    output_figure = output_entry + '_output_figure'
    output_tabular = output_entry + '_output_analysis.csv'
    output_full_table = output_entry + '_full_comparison_table.csv'
    
    is_gzip = gzip_flag or input_sequence.endswith('.gz') or output_entry.endswith('.gz')
    if output_entry.endswith('.gz'):
        output_entry = output_entry[:-3]
        
    if is_gzip:
        output_sequence = output_entry + '_output_clean.fasta.gz'
        output_rejected = output_entry + '_output_rejected.fasta.gz'
    else:
        output_sequence = output_entry + '_output_clean.fasta'
        output_rejected = output_entry + '_output_rejected.fasta'

    # Setup paths inside output_directory to avoid cluttering the workspace
    log_file_path = os.path.join(output_directory, output_entry + '_output.log')
    gap_plot_path = os.path.join(output_directory, output_entry + '_gap_plot.pdf')
    median_plot_path = os.path.join(output_directory, output_entry + '_median_plot.pdf')
    sampling_plot_path = os.path.join(output_directory, output_entry + '_sampling_trials_plot.pdf')

    # Logging and streaming to console
    mylogs = logging.getLogger(__name__)
    mylogs.setLevel(logging.INFO)
    
    # Avoid duplicate handlers if script is run multiple times in same session
    if not mylogs.handlers:
        stream = logging.StreamHandler()
        stream.setLevel(logging.INFO)
        streamformat = logging.Formatter("%(message)s")
        stream.setFormatter(streamformat)
        mylogs.addHandler(stream)

        file = logging.FileHandler(log_file_path)
        mylogs.addHandler(file)
    
    mylogs.info('\nSequenceBouncer: A method to remove outlier entries from a multiple sequence alignment\n')
    mylogs.info('Cory Dunn')
    mylogs.info('Contact : cory.david.dunn@gmail.com')
    mylogs.info('Version: ' + version)
    mylogs.info(f'Running with {cores} CPU core(s)/thread(s).')
    
    gpu_active = target_device is not None and device_flag in ('cuda', 'mps')
    if gpu_active:
        mylogs.info(f'GPU Acceleration Enabled ({target_device.type.upper()} available)')
    elif torch_available:
        mylogs.info('PyTorch CPU Acceleration Enabled')
    else:
        mylogs.info('Vectorized NumPy Mode (Standard)')
    mylogs.info('Please cite DOI: 10.1101/2020.11.24.395459')
    mylogs.info('___\n')
    
    start_time = time.time() 

    # Load sequences (support gzip)
    alignment_record_name_list = []
    if input_sequence.endswith('.gz'):
        import gzip
        with gzip.open(input_sequence, "rt") as f:
            records = list(SeqIO.parse(f, "fasta"))
    else:
        with open(input_sequence, "r") as f:
            records = list(SeqIO.parse(f, "fasta"))
            
    alignment_record_name_list = [record.name for record in records]
    depth_of_alignment = len(alignment_record_name_list)

    record_sequence_trial_results = pd.DataFrame(alignment_record_name_list, columns=['Accession'])

    length_of_alignment = len(records[0].seq)

    mylogs.info('Generating sequence dataframe.')
    
    # Vectorized fast sequence loading with progress bar
    sequence_dataframe = np.array([
        np.frombuffer(str(record.seq).lower().encode('ascii'), dtype=np.uint8) 
        for record in tqdm(records, desc="Loading sequences")
    ], dtype=np.uint8)

    # Release memory wrapper objects early
    del records
    gc.collect()

    # Initial properties
    unique_chars = np.unique(sequence_dataframe)
    N = depth_of_alignment

    mylogs.info('Calculating Shannon entropy values and gap metrics across all input sequences.')
    
    # Vectorized gap and shannon entropy calculation
    counts = np.zeros((len(unique_chars), length_of_alignment), dtype=np.int32)
    for idx, char in enumerate(unique_chars):
        counts[idx, :] = (sequence_dataframe == char).sum(axis=0)

    fractions = counts / depth_of_alignment
    nonzero_mask = fractions > 0
    log2_fractions = np.zeros_like(fractions)
    np.log2(fractions, where=nonzero_mask, out=log2_fractions)
    entropy_vector = -np.sum(fractions * log2_fractions, axis=0)

    if 45 in unique_chars:
        gap_idx = np.where(unique_chars == 45)[0][0]
        gap_record = fractions[gap_idx, :].tolist()
    else:
        gap_record = [0.0] * length_of_alignment

    # Plot gap fractions for alignment positions
    gap_record_sorted = sorted(gap_record)
    plotgaplistrange = np.arange(len(gap_record_sorted))
    plt.plot(plotgaplistrange, gap_record_sorted, 'o', ms=1, c='slategrey')
    plot_cutoff_label = 'Selected gap fraction cut-off: ' + str(gap_value_cutoff/100) + ' (' + str(gap_value_cutoff) + ' %)'
    plt.axhline(y=gap_value_cutoff/100, color='teal', linestyle='--', label=plot_cutoff_label)
    plt.xlabel('Input alignment column', fontsize=8)
    plt.ylabel('Gap fraction', fontsize=8)
    plt.legend()
    plt.savefig(gap_plot_path, format="pdf", bbox_inches="tight")
    mylogs.info('Printed gap distribution of input alignment to file: ' + gap_plot_path)
    plt.close()

    # Generate boolean based upon gap values
    gap_value_cutoff_float = float(gap_value_cutoff/100)
    gap_percent_bool_series_remove = np.array(gap_record) > gap_value_cutoff_float
    gap_percent_bool_index_remove = np.where(gap_percent_bool_series_remove)[0]

    # Remove gapped positions
    keep_indices = np.where(~gap_percent_bool_series_remove)[0]
    
    if len(keep_indices) == 0:
        mylogs.info('All columns were removed as gaps.')
        mylogs.info('Choose a larger value for --gap_percent_cut to continue.')
        exit()
    
    max_entropy_before_gaps = np.max(entropy_vector)
    mylogs.info('Maximum Shannon entropy alignment score before gap % considered: ' + str(round(max_entropy_before_gaps, 2)))
    
    entropy_gap_considered = entropy_vector[keep_indices]
    max_entropy_after_gaps = np.max(entropy_gap_considered)
    mylogs.info('Maximum Shannon entropy alignment score after gap % considered: ' + str(round(max_entropy_after_gaps, 2)))

    weights_vector = None
    if weighted:
        weights_vector = 1.0 / (entropy_gap_considered + entropy_smoothing)
        mylogs.info(f'Position-based weights initialized (Inverse Shannon entropy weighting enabled, smoothing={entropy_smoothing}).')

    mylogs.info('Removing gapped positions from analysis set.')
    sequence_dataframe_gap_considered = sequence_dataframe[:, keep_indices]
    
    mylogs.info("Elapsed time: ~ " + str(int(time.time() - start_time)) + " seconds.")
    mylogs.info('Alignment positions analyzed after ' + str(gap_value_cutoff) + '% gap cutoff: ' + str(len(keep_indices)))

    mylogs.info('Preparing sequences for comparison.')

    # -------- MEMORY CHECK AFTER GAP FILTERING --------
    unique_chars_filtered = np.unique(sequence_dataframe_gap_considered)
    K = len(unique_chars_filtered)
    L_filtered = len(keep_indices)
    
    gpu_mem, cpu_mem = get_available_memory()
    b_term = L_filtered * K
    
    # Calculate safe limits for GPU and CPU
    gpu_safe_N = 0
    if gpu_active and gpu_mem is not None:
        gpu_safe_N = int((-b_term + math.sqrt(b_term**2 + gpu_mem / 3.0)) / 2.0)
        gpu_safe_N = max(50, gpu_safe_N)
        
    cpu_safe_N = int((-b_term + math.sqrt(b_term**2 + cpu_mem / 3.0)) / 2.0)
    cpu_safe_N = max(50, cpu_safe_N)
    
    # Determine execution behavior based on mode, subsample_size, and available memory
    force_cpu_execution = (device_flag == 'cpu')
    
    # Check if GPU is target and check availability
    if device_flag in ('cuda', 'mps', 'gpu'):
        if not gpu_active:
            mylogs.error(f"\n[ERROR] Device '{device_flag}' is not available on this machine.")
            exit(1)

    # Initial mode determination
    mode = mode_flag
    if number_in_small_test > 0:
        if number_in_small_test == N:
            mode = 'full'
        else:
            mode = 'sample'

    # Determine device and feasibility for Full Analysis
    # We prefer GPU if available and it has enough memory for N
    # Otherwise CPU if CPU RAM is sufficient
    full_fits_on_gpu = gpu_active and (N <= gpu_safe_N)
    full_fits_on_cpu = N <= cpu_safe_N

    # Auto-resolution of mode and device
    chosen_mode = mode
    if mode == 'auto':
        if gpu_active and not force_cpu_execution:
            if full_fits_on_gpu:
                chosen_mode = 'full'
            elif full_fits_on_cpu:
                chosen_mode = 'full'
                force_cpu_execution = True
                mylogs.info(f"[INFO] Full analysis fits in CPU RAM ({cpu_mem / (1024**3):.2f} GB) but exceeds GPU VRAM. Running full analysis on CPU.")
            else:
                chosen_mode = 'sample'
                mylogs.info(f"[INFO] Full analysis exceeds system memory limit of N={max(gpu_safe_N, cpu_safe_N)}. Automatically running in sampling mode.")
        else:
            # CPU only
            if full_fits_on_cpu:
                chosen_mode = 'full'
            else:
                chosen_mode = 'sample'
                mylogs.info(f"[INFO] Full analysis exceeds CPU RAM limit of N={cpu_safe_N}. Automatically running in sampling mode.")

    # Now handle parameters based on the chosen mode
    if chosen_mode == 'full':
        # Warn if memory safety checks are violated (but allow override if user explicitly set mode='full' or number_in_small_test=N)
        if force_cpu_execution:
            if N > cpu_safe_N:
                if mode == 'full':
                    mylogs.warning(f"\n[WARNING] Forced full analysis requires N={N}, which exceeds estimated safe CPU memory limit of {cpu_safe_N}.")
                    mylogs.warning("          Execution will proceed on CPU, but may cause system slowdown or OOM crash.")
                else:
                    mylogs.error(f"\n[ERROR] Available CPU RAM ({cpu_mem / (1024**3):.2f} GB) is insufficient for a full analysis (requires N={N}).")
                    exit(1)
        else:
            # Running on GPU
            if N > gpu_safe_N:
                if mode == 'full':
                    mylogs.warning(f"\n[WARNING] Forced full analysis requires N={N}, which exceeds estimated safe GPU VRAM limit of {gpu_safe_N}.")
                    mylogs.warning("          Execution will proceed on GPU, but may crash with CUDA Out Of Memory.")
                elif N <= cpu_safe_N:
                    # Drop to CPU
                    force_cpu_execution = True
                    mylogs.warning(f"\n[WARNING] Available GPU memory ({gpu_mem / (1024**3):.2f} GB) is insufficient for a full analysis (requires N={N}).")
                    mylogs.warning(f"          Dropping down to CPU-based execution to perform the full analysis using available CPU RAM ({cpu_mem / (1024**3):.2f} GB).\n")
                else:
                    mylogs.error(f"\n[ERROR] Available memory is insufficient for a full analysis (requires N={N}).")
                    exit(1)

        number_in_small_test = N
        min_trials_for_each_sequence = 1
        is_full_analysis = True

    else:
        # sampling mode
        is_full_analysis = False
        active_safe_N = gpu_safe_N if (gpu_active and not force_cpu_execution) else cpu_safe_N
        
        # Calculate optimal subsample size if not explicitly set by the user
        if number_in_small_test <= 0:
            # Target sizes: GPU -> 2000, CPU -> 1000
            target_n = 2000 if (gpu_active and not force_cpu_execution) else 1000
            number_in_small_test = min(N, active_safe_N, target_n)
            # Ensure it is at least 50 (or N if N < 50)
            number_in_small_test = max(50 if N >= 50 else N, number_in_small_test)
            mylogs.info(f"[INFO] Automatically selected optimal subsample size: {number_in_small_test} (based on system resources).")
        else:
            # User overrode subsample size, check if it's safe
            if number_in_small_test > active_safe_N:
                if gpu_active and not force_cpu_execution and number_in_small_test <= cpu_safe_N:
                    if device_flag in ('cuda', 'mps', 'gpu'):
                        mylogs.warning(f"\n[WARNING] Requested subsample size {number_in_small_test} exceeds safe GPU limit of {gpu_safe_N} but --device '{device_flag}' is forced.")
                    else:
                        force_cpu_execution = True
                        mylogs.warning(f"\n[WARNING] Requested subsample size {number_in_small_test} exceeds GPU safe limit of {gpu_safe_N}.")
                        mylogs.warning(f"          Dropping down to CPU-based sampling execution using CPU RAM.\n")
                else:
                    mylogs.warning(f"\n[WARNING] Requested subsample size {number_in_small_test} exceeds safe limit of {active_safe_N}.")
                    mylogs.warning("          This may result in an Out-Of-Memory (OOM) crash.\n")

        # Calculate optimal trials if not explicitly set by the user
        if min_trials_for_each_sequence <= 0:
            # Heuristic for trials based on subsample size
            if number_in_small_test >= 1500:
                min_trials_for_each_sequence = 5
            elif number_in_small_test >= 1000:
                min_trials_for_each_sequence = 6
            elif number_in_small_test >= 500:
                min_trials_for_each_sequence = 8
            elif number_in_small_test >= 200:
                min_trials_for_each_sequence = 10
            elif number_in_small_test >= 100:
                min_trials_for_each_sequence = 12
            else:
                min_trials_for_each_sequence = 15
            mylogs.info(f"[INFO] Automatically selected optimal trials count: {min_trials_for_each_sequence} (for subsample size {number_in_small_test}).")

    mylogs.info("\nAnalyzing '" + input_sequence + "'.")
    mylogs.info('Flags are --IQR_coefficient: ' + str(multiplier_on_interquartile_range) + ', -subsample_size: ' + str(number_in_small_test) + ', --gap_percent_cut: ' + str(gap_value_cutoff))
    if min_trials_for_each_sequence != 1:
        mylogs.info('          --stringency: ' + str(stringency_flag) + ', --trials: ' + str(min_trials_for_each_sequence) + ', --random_seed: ' + str(seed))

    mylogs.info('Input alignment length is: ' + str(length_of_alignment) + ' characters.')
    mylogs.info("Input alignment depth is: " + str(depth_of_alignment) + " sequences.")

    # Convert to pandas-like structure for metadata tracking
    original_row_indices = np.arange(depth_of_alignment)

    # Clear out unused items from memory
    del sequence_dataframe
    gc.collect()

    record_sequence_trial_results['Total_trials'] = 0
    record_sequence_trial_results['Outlier_instances'] = 0

    mylogs.info('Beginning sequence trials.')
    
    # Initialize random state generator
    rng = np.random.default_rng(seed)

    if depth_of_alignment // number_in_small_test == depth_of_alignment / number_in_small_test:
        times_to_sample_max_keep = (depth_of_alignment // number_in_small_test)
    else:
        times_to_sample_max_keep = (depth_of_alignment // number_in_small_test) + 1

    for trial in range(min_trials_for_each_sequence):
        if min_trials_for_each_sequence > 1:
            mylogs.info("Trial: " + str(trial+1) + " of " + str(min_trials_for_each_sequence))
        
        # Shuffle row indices reproducibly using rng
        shuffled_rows = rng.permutation(depth_of_alignment)
        shuffled_seq_matrix = sequence_dataframe_gap_considered[shuffled_rows]
        shuffled_orig_indices = original_row_indices[shuffled_rows]

        # Setup sample range with tqdm progress bar
        sample_range = range(times_to_sample_max_keep)
        if number_in_small_test != depth_of_alignment:
            sample_range = tqdm(sample_range, desc=f"Trial {trial+1}/{min_trials_for_each_sequence}", leave=False)

        for j in sample_range: 
            # Find which sequences need to be tested (Total_trials < trial + 1)
            current_trials = record_sequence_trial_results.loc[shuffled_orig_indices, 'Total_trials'].values
            eligible_mask = current_trials <= trial
            
            eligible_seq = shuffled_seq_matrix[eligible_mask]
            eligible_indices = shuffled_orig_indices[eligible_mask]
            
            if len(eligible_seq) == 0:
                continue
                
            number_to_choose = min(number_in_small_test, len(eligible_seq))
            
            table_sample_numpy = eligible_seq[:number_to_choose]
            table_sample_indices = eligible_indices[:number_to_choose]
            
            # Distance computation using optimized function
            transformed_distances = compute_distances(table_sample_numpy, max_entropy_after_gaps, weights=weights_vector, force_cpu=force_cpu_execution)
            
            # Median along axis 1 (row-wise)
            medians = np.median(transformed_distances, axis=1)
            
            # IQR calculations
            q25, q75 = np.percentile(medians, [25, 75])
            iqr = q75 - q25
            CIQR = iqr * multiplier_on_interquartile_range
            upper_cutoff = q75 + CIQR
            
            # Update trials
            record_sequence_trial_results.loc[table_sample_indices, 'Total_trials'] += 1
            
            # Mark outliers
            outlier_mask = medians > upper_cutoff
            outlier_indices = table_sample_indices[outlier_mask]
            record_sequence_trial_results.loc[outlier_indices, 'Outlier_instances'] += 1

            # If full matrix analysis and first trial, generate outputs
            if number_in_small_test == depth_of_alignment:
                # Plot comparisons
                medians_sorted = sorted(medians.tolist())
                plotlistrange = np.arange(len(medians_sorted))
                plt.plot(plotlistrange, medians_sorted, 'o', ms=1, c='darkgreen')
                plot_cutoff_label = 'Selected IQR cut-off:  ' + str(multiplier_on_interquartile_range)
                plt.axhline(y=upper_cutoff, color='red', linestyle='--', label=plot_cutoff_label)
                plt.xlabel('Sequence', fontsize=8)
                plt.ylabel('Median across pairwise comparisons', fontsize=8)
                plt.legend()
                plt.savefig(median_plot_path, format="pdf", bbox_inches="tight")
                mylogs.info('Printed median values of sequence comparisons from full analysis to file ' + median_plot_path)
                plt.close()

                # Generate full comparison table
                mylogs.info('Cut-off value for median taken across comparisons (full-alignment pairwise analysis): ' + str(round(upper_cutoff,1)))
                
                # Make DataFrame with accessions
                entropy_DF = pd.DataFrame(transformed_distances, index=alignment_record_name_list, columns=alignment_record_name_list)
                entropy_DF.insert(0, 'Median_across_pairwise_comparisons', medians)
                entropy_DF.to_csv(os.path.join(output_directory, output_full_table))
        
        mylogs.info("Elapsed time: ~ " + str(int(time.time() - start_time)) + " seconds.")
        mylogs.info("Estimated total time for analysis: ~ " + str(int(((time.time() - start_time))/(1+trial)*min_trials_for_each_sequence)) + " seconds.")

    # Free memory cache from GPU
    if torch_available:
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            torch.mps.empty_cache()

    # Process and write output FASTA files
    record_sequence_trial_results['Fraction_positive'] = record_sequence_trial_results['Outlier_instances'] / record_sequence_trial_results['Total_trials']
    
    # Reload sequences from FASTA (support gzip)
    if input_sequence.endswith('.gz'):
        import gzip
        with gzip.open(input_sequence, "rt") as f:
            records = list(SeqIO.parse(f, "fasta"))
    else:
        with open(input_sequence, "r") as f:
            records = list(SeqIO.parse(f, "fasta"))
    record_seqs = [str(record.seq) for record in records]
    frame = { 'Accession': alignment_record_name_list, 'Sequence': record_seqs }
    FASTA_output_unclean_DF = pd.DataFrame(frame)

    # Plot trials for sampling-based approach
    if number_in_small_test != depth_of_alignment:
        record_sequence_trial_results_list = record_sequence_trial_results['Fraction_positive'].tolist()
        record_sequence_trial_results_list.sort()
        plottrialrange = np.arange(len(record_sequence_trial_results_list))
        plt.plot(plottrialrange, record_sequence_trial_results_list, 'o', ms=1, c='deepskyblue')
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['Arial']
        plt.xlabel('Sequence', fontsize=8)
        plt.ylabel('Fraction of times sequence called aberrant (In order of increasing positive calls)', fontsize=8)
        plt.savefig(sampling_plot_path, format="pdf", bbox_inches="tight")
        mylogs.info('Printed plot of sampling trial results to file: ' + sampling_plot_path)
        plt.close()

    # Apply stringency settings
    if stringency_flag == 1:
        clean_mask = record_sequence_trial_results['Fraction_positive'] != 1
    elif stringency_flag == 2:
        clean_mask = record_sequence_trial_results['Fraction_positive'] <= 0.5
    else: # stringency 3
        clean_mask = record_sequence_trial_results['Fraction_positive'] == 0

    FASTA_output_clean_DF = FASTA_output_unclean_DF.loc[clean_mask]
    FASTA_output_reject_DF = FASTA_output_unclean_DF.loc[~clean_mask]

    # Save clean FASTA file
    mylogs.info('Writing cleaned alignment as FASTA.')
    mylogs.info(FASTA_output_clean_DF)
    FASTA_output_final_acc_list = FASTA_output_clean_DF['Accession'].tolist()
    FASTA_output_final_seq_list = FASTA_output_clean_DF['Sequence'].tolist()

    def open_fasta_out(filepath):
        if filepath.endswith('.gz'):
            import gzip
            return gzip.open(filepath, "wt")
        return open(filepath, "w")

    with open_fasta_out(os.path.join(output_directory, output_sequence)) as ofile:
        for acc, seq in zip(FASTA_output_final_acc_list, FASTA_output_final_seq_list):
            ofile.write(f">{acc}\n{seq}\n")

    # Save rejected FASTA file
    mylogs.info('Writing rejected sequences to FASTA.')
    mylogs.info(FASTA_output_reject_DF)
    FASTA_output_rejected_acc_list = FASTA_output_reject_DF['Accession'].tolist()
    FASTA_output_rejected_seq_list = FASTA_output_reject_DF['Sequence'].tolist()
        
    with open_fasta_out(os.path.join(output_directory, output_rejected)) as ofile:
        for acc, seq in zip(FASTA_output_rejected_acc_list, FASTA_output_rejected_seq_list):
            ofile.write(f">{acc}\n{seq}\n")

    # Save analysis file as CSV
    record_sequence_trial_results['Selected_for_retention'] = np.where(clean_mask, True, False)

    if min_trials_for_each_sequence == 1:
        record_sequence_trial_results = record_sequence_trial_results.drop(columns=['Total_trials', 'Outlier_instances', 'Fraction_positive'])

    record_sequence_trial_results.to_csv(os.path.join(output_directory, output_tabular), index=False)

    mylogs.info("Analysis complete.")
    mylogs.info("Total time for analysis: ~ " + str(int(((time.time() - start_time))/(1+trial)*min_trials_for_each_sequence)) + " seconds.")
