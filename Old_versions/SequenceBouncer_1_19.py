#!/usr/bin/env python
# coding: utf-8 

# Funding received from the Sigrid Jusélius Foundation contributed to the development of this software.
# Author: Cory Dunn
# Institution: University of Helsinki
# Author Email: cory.dunn@helsinki.fi
# Version: 1.19
version = '1.19'
# License: GPLv3

from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import SeqIO
import pandas as pd
import numpy as np
import time
import sys
import math
import argparse
import gc
import random

print('\nSequenceBouncer: A method to remove outlier entries from a multiple sequence alignment\n')
print('Cory Dunn')
print('University of Helsinki')
print('cory.dunn@helsinki.fi')
print('Version: ' + version)
print('Please cite DOI: 10.1101/2020.11.24.395459')
print('___\n')

# Load files, receive parameters, and provide assistance

ap = argparse.ArgumentParser()
ap.add_argument('-i','--input_file',required=True,type=str,help='Input file in FASTA format.\n')
ap.add_argument('-o','--output_file',required=False,type=str,default='X',help="Output filename [do not include extensions] (default will be 'input_file.ext').\n")
ap.add_argument('-g','--gap_percent_cut',required=False,type=float,default=2.0,help='For columns with a greater fraction of gaps than the selected value, expressed in percent, data will be ignored in calculations (default is 2).\n')
ap.add_argument('-k','--IQR_coefficient',required=False,type=float,default=1.0,help='Coefficient multiplied by the interquartile range that helps to define an outlier sequence (default is 1.0).\n')
ap.add_argument('-n','--subsample_size',required=False,type=int,default=0,help='|> Available for large alignments | The size of a single sample taken from the full dataset (default is entire alignment, but try a subsample size of 50 or 100 for large alignments).\n')
ap.add_argument('-t','--trials',required=False,type=int,default=1,help='|> Available for large alignments | Number of times each sequence is sampled and tested (default is to examine all sequences in one single trial, but 5 or 10 trials may work well when subsamples are taken from large alignments).\n')
ap.add_argument('-s','--stringency',required=False,type=int,default=2,help='|> Available for large alignments | 1: Minimal stringency 2: Moderate stringency 3: Maximum stringency (default is moderate stringency).\n')
ap.add_argument('-r','--random_seed',required=False,type=int,default=random.randint(0,1000),help='Random seed (integer) to be used during a sampling-based approach (default is that the seed is randomly selected). The user can use this seed to obtain reproducible output and should note it in their publications. \n')

args = vars(ap.parse_args())
input_sequence = args['input_file']
stringency_flag = args['stringency']
min_trials_for_each_sequence = args['trials']
multiplier_on_interquartile_range = args['IQR_coefficient']
number_in_small_test = args['subsample_size']
gap_value_cutoff = args['gap_percent_cut']
output_entry = args['output_file']
seed = args['random_seed']

# Prepare output filenames

if output_entry == 'X':
    sep = '.'
    input_strip = input_sequence.split(sep, 1)[0]
    output_figure = input_strip + '_output_figure'
    output_sequence = input_strip + '_output_clean.fasta'
    output_rejected = input_strip + '_output_rejected.fasta'
    output_tabular = input_strip + '_output_analysis.csv'
    output_full_table = input_strip + '_full_comparison_table.csv'
elif output_entry != 'X':
    output_figure = output_entry + '_output_figure'
    output_sequence = output_entry + '_output_clean.fasta'
    output_rejected = output_entry + '_output_rejected.fasta'
    output_tabular = output_entry + '_output_analysis.csv'
    output_full_table = output_entry + '_full_comparison_table.csv'

# Start timer

start_time = time.time() 

# Initialize

alignment_record_name_list = []

for record in SeqIO.parse(input_sequence,"fasta"):
        alignment_record_name_list.append(record.name)

depth_of_alignment = (len(alignment_record_name_list))

record_sequence_trial_results = pd.DataFrame(alignment_record_name_list, columns=['Accession'])

if number_in_small_test == 0:
    number_in_small_test = depth_of_alignment

if number_in_small_test == depth_of_alignment:
    min_trials_for_each_sequence = 1

print("Analyzing '" + input_sequence + "'.")
print('Flags are --IQR_coefficient: ' + str(multiplier_on_interquartile_range) + ', -subsample_size: ' + str(number_in_small_test) + ', --gap_percent_cut: ' + str(gap_value_cutoff))
if min_trials_for_each_sequence != 1:
    print('          --stringency: ' + str(stringency_flag) + ', --trials: ' + str(min_trials_for_each_sequence) + ', --random_seed: ' + str(seed))

length_of_alignment = len(list(record.seq))

print('Input alignment length is: ' + str(length_of_alignment) + ' characters.')
print("Input alignment depth is: " + str(depth_of_alignment) + " sequences.")

# Load sequences from alignment into list and control case

print('Generating sequence dataframe.')

record_x_toward_seq_dataframe = []
sequence_records = []

for record_x in SeqIO.parse(input_sequence,"fasta"):
    record_x_toward_seq_dataframe = list(record_x.seq)
    record_x_toward_seq_dataframe_lower = [x.lower() for x in record_x_toward_seq_dataframe] 
    record_x_toward_seq_dataframe_ASCII = [ord(x) for x in record_x_toward_seq_dataframe_lower]
    sequence_records.append(record_x_toward_seq_dataframe_ASCII)

# Generate dataframe of alignment from list

sequence_dataframe = pd.DataFrame(sequence_records)
sequence_dataframe = sequence_dataframe.astype('int8')

# Calculate Shannon entropy and fraction of column gapped

print('Calculating Shannon entropy values and gap metrics across all input sequences.')
entropy_record = []
gap_record = []
sequence_columns = len(sequence_dataframe.axes[1])
for i in range(sequence_columns):
    column_fractions_S = sequence_dataframe[i].value_counts(normalize='True')
    shannon_entropy_column = 0
    gap_fraction_column = 0
    for character, val in  column_fractions_S.iteritems():
        shannon_entropy_column +=  val * math.log(val,2)
        if character == 45:
            gap_fraction_column = val
    shannon_entropy_column *= -1
    entropy_record.append(shannon_entropy_column)
    gap_record.append(gap_fraction_column)
entropylist_S = pd.Series(entropy_record)
gap_fraction_S = pd.Series(gap_record)

# Generate boolean based upon gap values

gap_value_cutoff_float = float(gap_value_cutoff/100)
gap_percent_bool_series_remove = gap_fraction_S > gap_value_cutoff_float
gap_percent_bool_index_remove = gap_percent_bool_series_remove[gap_percent_bool_series_remove].index

# Remove gapped positions

entropylist_S_gap_considered = entropylist_S.drop(gap_percent_bool_index_remove)
max_entropy_before_gaps = pd.Series.max(entropylist_S)
print('Maximum Shannon entropy alignment score before gap % considered: ', f'{max_entropy_before_gaps: .2f}')
max_entropy_after_gaps = pd.Series.max(entropylist_S_gap_considered)
print('Maximum Shannon entropy alignment score after gap % considered: ', f'{max_entropy_after_gaps: .2f}')

entropy_record_numpy = entropylist_S_gap_considered.to_numpy()
entropy_record_numpy.shape = (-1,len(entropylist_S_gap_considered))

print('Removing gapped positions from analysis set.')
sequence_dataframe_gap_considered = sequence_dataframe.drop(gap_percent_bool_index_remove,axis=1)
print("Elapsed time: ~ " + str(int(time.time() - start_time)) + " seconds.")
print('Alignment positions analyzed after ' + str(gap_value_cutoff) + '% gap cutoff: ' + str(length_of_alignment-len(gap_percent_bool_index_remove)))

print('Preparing sequences for comparison.')

# Comparison time warning

comparison_time_full_table_seconds = depth_of_alignment * depth_of_alignment * (length_of_alignment-len(gap_percent_bool_index_remove)) * 3.14E-8
if comparison_time_full_table_seconds > 1800 and number_in_small_test == depth_of_alignment:
    print('\n***WARNING: An input alignment of this size may take a considerable amount of time')
    print('   if all pairwise sequence comparisons are performed.')
    print('   A sampling-based approach may be considered.')
    print('   For a sampling-based approach, take advantage of the -n, -t, and -s flags.\n')

# Clear out unused items from memory

del sequence_records
del sequence_dataframe
gc.collect()

# Prepare dataframe for storage of trial results (these columns are stripped away later if only one trial is performed)

record_sequence_trial_results['Total_trials'] = 0
record_sequence_trial_results['Outlier_instances'] = 0

# Set trial counter

print('Beginning sequence trials.')
trial_count = 0

# Avoid empty source dataframe

if depth_of_alignment//number_in_small_test == depth_of_alignment/number_in_small_test:
    times_to_sample_max_keep = (depth_of_alignment//number_in_small_test)
else:
    times_to_sample_max_keep = (depth_of_alignment//number_in_small_test) + 1

# Define the calculation engine

def engine():
    for counter_x in range(table_sample_numpy_rows):
                counter_x_numpy_row = table_sample_numpy[counter_x:(counter_x+1),:]
                if depth_of_alignment < 1000 and ((counter_x+1)/25) == ((counter_x+1)//25):
                    print('\rSequences analyzed: '+str(counter_x+1))
                elif depth_of_alignment < 10000 and ((counter_x+1)/250) == ((counter_x+1)//250):
                    print('\rSequences analyzed: '+str(counter_x+1))
                elif depth_of_alignment < 100000 and ((counter_x+1)/2500) == ((counter_x+1)//2500):
                    print('\rSequences analyzed: '+str(counter_x+1))
                for counter_y in range((counter_x+1)):
                    counter_y_numpy_row = table_sample_numpy[counter_y:(counter_y+1),:]
                    comparison_bool_series_match = counter_x_numpy_row == counter_y_numpy_row
                    comparison_bool_series_NOT_match = counter_x_numpy_row != counter_y_numpy_row
                    entropy_record_match = entropy_record_numpy[(comparison_bool_series_match)]
                    entropy_record_NOT_match = entropy_record_numpy[(comparison_bool_series_NOT_match)]
                    match_entropy_total = entropy_record_match.sum(axis=0)
                    NOT_match_entropy_minus_max_entropy = entropy_record_NOT_match - max_entropy_after_gaps
                    NOT_match_entropy_total =  NOT_match_entropy_minus_max_entropy.sum(axis=0)
                    total_entropy_recorded = match_entropy_total + NOT_match_entropy_total
                    entropy_array[counter_x, counter_y] = total_entropy_recorded
                    entropy_array[counter_y, counter_x] = total_entropy_recorded
    return entropy_array

for trial in range(min_trials_for_each_sequence):
    
    if min_trials_for_each_sequence > 1:
        print("Trial: " + str(trial+1) + " of " + str(min_trials_for_each_sequence))
    
    sequence_dataframe_gap_considered = sequence_dataframe_gap_considered.sample(frac=1,random_state = seed) # shuffle master sequence dataframe, use user-defined or standard random seed
    sequence_dataframe_gap_considered_max_keep = sequence_dataframe_gap_considered # copy shuffled version for work below

    for j in range(times_to_sample_max_keep): 
        if number_in_small_test != depth_of_alignment and (j+1)//50 == (j+1)/50:
            print('\rSample: '+str((j+1)) + ' of ' +str(times_to_sample_max_keep) + ' | Trial: ' + str(trial+1))
        max_times_tested = record_sequence_trial_results['Total_trials'].max()
        
        if max_times_tested > trial:
            
            sequence_dataframe_gap_considered_max_keep = sequence_dataframe_gap_considered.loc[record_sequence_trial_results['Total_trials'] != max_times_tested]
    
        length_sequence_dataframe_gap_considered_max_keep = len(sequence_dataframe_gap_considered_max_keep)
        
        if length_sequence_dataframe_gap_considered_max_keep >= number_in_small_test:
            number_to_choose = number_in_small_test
    
        elif length_sequence_dataframe_gap_considered_max_keep < number_in_small_test:
            number_to_choose = length_sequence_dataframe_gap_considered_max_keep
        
        table_sample = sequence_dataframe_gap_considered_max_keep.iloc[0:number_to_choose,:]
        table_sample_numpy = table_sample.to_numpy() # convert pandas dataframe to numpy array
        table_sample_numpy = table_sample_numpy.astype(np.int8) # change datatype in an attempt to reduce memory imprint
        table_sample_numpy_rows, table_sample_numpy_columns = table_sample_numpy.shape
    
    # Initiate numpy array for entropy calculation values
    
        entropy_array = np.empty((number_to_choose,number_to_choose),dtype=float)
        entropy_array[:] = np.nan
                          
    # Calculations of match or not, and sum entropy values

        engine()

        entropy_DF = pd.DataFrame(entropy_array,index=np.arange(number_to_choose), columns=np.arange(number_to_choose))
        maximum_entropy_sum = entropy_DF.max(skipna=True)
        entropy_DF -= maximum_entropy_sum
        entropy_DF *= -1.0
        entropy_DF.columns = table_sample.index
        entropy_DF.index = table_sample.index
        entropy_DF_analysis_empty = np.empty((number_to_choose,1))
        entropy_DF_analysis_empty[:] = np.nan
        entropy_DF_analysis = pd.DataFrame(data = entropy_DF_analysis_empty, index=entropy_DF.index, columns=['Median'])
        
        for z in entropy_DF:
            entropy_DF_analysis.loc[z,'Median'] = entropy_DF.loc[z,:].median(skipna=True)
       
        record_sequence_trial_results.loc[entropy_DF_analysis.index,'Total_trials'] += 1

# Calculate interquartile range and outlier cutoff

        entropy_DF_analysis_values_list = entropy_DF_analysis.values.tolist()
        q25, q75 = np.nanpercentile(entropy_DF_analysis_values_list, 25), np.nanpercentile(entropy_DF_analysis_values_list, 75)
        iqr = q75 - q25

        CIQR = iqr * multiplier_on_interquartile_range
        lower_cutoff, upper_cutoff = q25 - CIQR, q75 + CIQR

# Identify the outlier sequences using the interquartile range cutoff

        entropy_DF_analysis_above_cutoff = entropy_DF_analysis > upper_cutoff
        entropy_median_too_high = entropy_DF_analysis_above_cutoff.loc[entropy_DF_analysis_above_cutoff['Median'] == True]
        record_sequence_trial_results.loc[entropy_median_too_high.index,'Outlier_instances'] += 1
    
    print("Elapsed time: ~ " + str(int(time.time() - start_time)) + " seconds.")
    print("Estimated total time for analysis: ~ " + str(int(((time.time() - start_time))/(1+trial)*min_trials_for_each_sequence)) + " seconds.")

# Print full distance matrix for analysis and generate a plot only if a single test of all sequences versus all sequences was performed

if number_in_small_test == depth_of_alignment:

    print('Cut-off value for median taken across comparisons (full-alignment pairwise analysis): ', f'{upper_cutoff: .1f}')
    entropy_DF.sort_index(axis=0,inplace=True,ascending=True)
    entropy_DF.sort_index(axis=1,inplace=True,ascending=True)
    entropy_DF.index = alignment_record_name_list
    entropy_DF.columns = alignment_record_name_list
    entropy_DF['Median_across_pairwise_comparisons'] = entropy_DF.median(axis=1) # add column calculating median across pairwise comparisons
    first_column = entropy_DF.pop('Median_across_pairwise_comparisons')
    entropy_DF.insert(0,'Median_across_pairwise_comparisons',first_column)
    entropy_DF.to_csv(output_full_table)
    
    # Print figure describing trials (NOT ACTIVE IN CURRENT VERSION)
    # import matplotlib.pyplot as plt
    # from matplotlib import rcParams
    # rcParams['font.family'] = 'sans-serif'
    # rcParams['font.sans-serif'] = ['Arial']
    # fig, ax = plt.subplots()
    # plt.figure(figsize=(8,8))
    # plt.xlabel('Sequence number within input file', fontsize=12)
    # plt.ylabel('Sequence number within input file', fontsize=12)
    # if depth_of_alignment < 2000:
    #    plotlabels = [i+1 for i in range(0,(depth_of_alignment-1),50)]
    # if depth_of_alignment >= 2000:
    #    plotlabels = [i+1 for i in range(0,(depth_of_alignment-1),500)]
    # plt.xticks(ticks=plotlabels,rotation=60)
    # plt.yticks(ticks=plotlabels)
    # plt.imshow(entropy_DF,cmap='cividis',interpolation='nearest')
    # plt.show
    # plt.colorbar()
    # plt.savefig(output_figure + '.pdf', dpi=600)

# Prepare dataframe to generate FASTA files

record_sequence_trial_results['Fraction_positive'] = record_sequence_trial_results['Outlier_instances'] / record_sequence_trial_results['Total_trials']
record_seq_convert_to_string = []

for record in SeqIO.parse(input_sequence,"fasta"):
    record_seq_convert_to_string.append(str(record.seq))
    
acc_records_S = pd.Series(alignment_record_name_list)
sequence_records_S = pd.Series(record_seq_convert_to_string)
    
frame = { 'Accession': acc_records_S, 'Sequence': sequence_records_S }
FASTA_output_unclean_DF = pd.DataFrame(frame) 

# Generating clean dataframes

if stringency_flag == 1: # minimal stringency
    FASTA_output_clean_DF = FASTA_output_unclean_DF.loc[record_sequence_trial_results['Fraction_positive'] != 1]
if stringency_flag == 2: # moderate stringency
    FASTA_output_clean_DF = FASTA_output_unclean_DF.loc[record_sequence_trial_results['Fraction_positive'] <= 0.5]
if stringency_flag == 3: # maximal stringency
    FASTA_output_clean_DF = FASTA_output_unclean_DF.loc[record_sequence_trial_results['Fraction_positive'] == 0]

# Generating rejection dataframes

if stringency_flag == 1: # minimal stringency
    FASTA_output_reject_DF = FASTA_output_unclean_DF.loc[record_sequence_trial_results['Fraction_positive'] == 1]
if stringency_flag == 2: # moderate stringency
    FASTA_output_reject_DF = FASTA_output_unclean_DF.loc[record_sequence_trial_results['Fraction_positive'] > 0.5]
if stringency_flag == 3: # maximal stringency
    FASTA_output_reject_DF = FASTA_output_unclean_DF.loc[record_sequence_trial_results['Fraction_positive'] != 0]

# Save clean FASTA file

print('Writing cleaned alignment as FASTA.')
print(FASTA_output_clean_DF)
FASTA_output_final_acc_list = FASTA_output_clean_DF.loc[:,'Accession'].values.tolist()
FASTA_output_final_seq_list = FASTA_output_clean_DF.loc[:,'Sequence'].values.tolist()

ofile = open(output_sequence, "w")
    
for seqi in range(len(FASTA_output_final_acc_list)):
    ofile.write(">" + FASTA_output_final_acc_list[seqi] + "\n" + FASTA_output_final_seq_list[seqi] + "\n")
ofile.close()

# Save rejected FASTA file

print('Writing rejected sequences to FASTA.')
print(FASTA_output_reject_DF)
FASTA_output_rejected_acc_list = FASTA_output_reject_DF.loc[:,'Accession'].values.tolist()
FASTA_output_rejected_seq_list = FASTA_output_reject_DF.loc[:,'Sequence'].values.tolist()
    
ofile = open(output_rejected, "w")
for seqi in range(len(FASTA_output_rejected_acc_list)):
    ofile.write(">" + FASTA_output_rejected_acc_list[seqi] + "\n" + FASTA_output_rejected_seq_list[seqi] + "\n")
ofile.close()

# Save analysis file as CSV

record_sequence_trial_results['Selected_for_retention'] = ''
if stringency_flag == 1: # minimal stringency
    record_sequence_trial_results['Selected_for_retention'] = np.where(record_sequence_trial_results['Fraction_positive'] != 1 ,True,False)

if stringency_flag == 2: # moderate stringency
    record_sequence_trial_results['Selected_for_retention'] = np.where(record_sequence_trial_results['Fraction_positive'] <= 0.5,True,False)

if stringency_flag == 3: # maximal stringency
    record_sequence_trial_results['Selected_for_retention'] = np.where(record_sequence_trial_results['Fraction_positive'] == 0 ,True,False)

if min_trials_for_each_sequence == 1:
    del record_sequence_trial_results['Total_trials']
    del record_sequence_trial_results['Outlier_instances']
    del record_sequence_trial_results['Fraction_positive']

record_sequence_trial_results.to_csv(output_tabular,index=False)

# Provide total runtime

print("Analysis complete.")
print("Total time for analysis: ~ " + str(int(((time.time() - start_time))/(1+trial)*min_trials_for_each_sequence)) + " seconds.")

# %%