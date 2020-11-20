#!/usr/bin/env python
# coding: utf-8

# Funding received from the Sigrid Jusélius Foundation contributed to the development of this software.
# Author: Cory Dunn
# Institution: University of Helsinki
# Author Email: cory.dunn@helsinki.fi
# Version: 1.17-A
version = '1.17-A'
# License: GPLv3

from Bio import AlignIO #import AlignIO package
from Bio.Align import AlignInfo #import AlignInfo
from Bio import SeqIO #import SeqIO package
import pandas as pd #import pandas
import numpy as np #import numpy
import time #import time
import sys
import argparse
#import matplotlib.pyplot as plt
#from matplotlib import rcParams
#rcParams['font.family'] = 'sans-serif'
#rcParams['font.sans-serif'] = ['Arial']
import gc

print('\nSequenceBouncer: A method to remove outlier entries from a multiple sequence alignment\n')
print('Cory Dunn')
print('University of Helsinki')
print('cory.dunn@helsinki.fi')
print('Version: ' + version +'\n')

# Load files, receive parameters, and provide assistance

ap = argparse.ArgumentParser()
ap.add_argument('-i','--input_file',required=True,type=str,help='Input file in FASTA format.\n')
ap.add_argument('-o','--output_file',required=False,type=str,default='X',help="Output filename [do not include extensions] (Default will be 'input_file.ext').\n")
ap.add_argument('-g','--gap_percent_cut',required=False,type=float,default=2.0,help='For columns with a greater fraction of gaps than the selected value (expressed in percent), data will be ignored in calculations (Default is 2).\n')
ap.add_argument('-k','--IQR_coefficient',required=False,type=float,default=1.0,help='Coefficient multiplied by the interquartile range that helps to define an outlier sequence (Default is 1.0).\n')
ap.add_argument('-n','--subsample_size',required=False,type=int,default=0,help='|> Available for large alignments | The size of a single sample taken from the full dataset (Default is entire alignment, but try a subsample size of 50 or 100 for large alignments).\n')
ap.add_argument('-t','--trials',required=False,type=int,default=1,help='|> Available for large alignments | Number of times each sequence is sampled and tested (Default is to examine all sequences in one single trial, but 5 or 10 trials may work well when subsamples are taken from large alignments).\n')
ap.add_argument('-s','--stringency',required=False,type=int,default=2,help='|> Available for large alignments | 1: Minimal stringency. 2: Moderate stringency. 3: Maximum stringency. (Default is moderate stringency).\n')

args = vars(ap.parse_args())
input_sequence = args['input_file']
stringency_flag = args['stringency']
min_trials_for_each_sequence = args['trials']
multiplier_on_interquartile_range = args['IQR_coefficient']
number_in_small_test = args['subsample_size']
gap_value_cutoff = args['gap_percent_cut']
output_entry = args['output_file']

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

# start timer

start_time = time.time() 

# initialize numpy array for sum entropy values

alignment_record_name_list = []
entropy_array = []

for record in SeqIO.parse(input_sequence,"fasta"):
        alignment_record_name_list.append(record.name)

depth_of_alignment = (len(alignment_record_name_list))

entropy_array = np.empty((depth_of_alignment,depth_of_alignment))
entropy_array[:] = np.nan

record_sequence_trial_results = pd.DataFrame(alignment_record_name_list, columns=['Accession'])

if number_in_small_test == 0:
    number_in_small_test = depth_of_alignment

if number_in_small_test == depth_of_alignment:
    min_trials_for_each_sequence = 1

print("Analyzing '" + input_sequence + "'.")
print('Flags are -k: ' + str(multiplier_on_interquartile_range) + ', -n: ' + str(number_in_small_test) + ', -g: ' + str(gap_value_cutoff))
if min_trials_for_each_sequence != 1:
    print('          -s: ' + str(stringency_flag) + ', -t: ' + str(min_trials_for_each_sequence))

# initialize numpy array for gap counts

length_of_alignment = len(list(record.seq))

gap_array = np.zeros((depth_of_alignment,length_of_alignment))

print('Alignment length is: ' + str(int(gap_array.size/depth_of_alignment)) + ' characters.')
print("Alignment depth is: " + str(depth_of_alignment) + " sequences.")

# test gap percentage at different alignment positions

print('Calculating gap metrics.')
count_a = -1

for record in SeqIO.parse(input_sequence,"fasta"):
    count_a += 1  
    alignment_record_sequence_list_for_gap_breakdown = list(record.seq)
    for idx, item in enumerate(alignment_record_sequence_list_for_gap_breakdown): 
        if alignment_record_sequence_list_for_gap_breakdown[idx] == '-':
            gap_array[count_a,idx] = 1
        else:
            gap_array[count_a,idx] = 0

del alignment_record_sequence_list_for_gap_breakdown

gap_array_column_sums = gap_array.sum(axis=0)

gap_array_column_sum_S = pd.Series(gap_array_column_sums)

# generate boolean based upon gap values

gap_percent = gap_array_column_sum_S/depth_of_alignment

gap_value_cutoff_float = float(gap_value_cutoff/100)

gap_percent_bool_series_remove = gap_percent > gap_value_cutoff_float

gap_percent_bool_index_remove = gap_percent_bool_series_remove[gap_percent_bool_series_remove].index

print('Positions analyzed after ' + str(gap_value_cutoff) + '% gap cutoff: ' + str(length_of_alignment-len(gap_percent_bool_index_remove)))

comparison_time_full_table_seconds = depth_of_alignment * depth_of_alignment * (length_of_alignment-len(gap_percent_bool_index_remove)) * 3.14E-8
if comparison_time_full_table_seconds > 1800 and number_in_small_test == depth_of_alignment:
    print('\n***WARNING: An input alignment of this size may take a considerable amount of time')
    print('   if all pairwise sequence comparisons are performed.')
    print('   A sampling-based approach may be considered.')
    print('   For a sampling-based approach, take advantage of the -n, -t, and -s flags.\n')

print('Calculating Shannon entropy values across all input sequences.')

# Calculate Shannon entropy from a MSA. Routine derived from ShannonMSA 1.0.0 by Joe R. J. Healey, University of Warwick
# under the GPLv3 license. Gaps and N's are included in the calculation.

msa = input_sequence
alnformat = "fasta"
verbose = 1

def parseMSA(msa, alnformat, verbose):
    from Bio import AlignIO
    alignment = AlignIO.read(msa, alnformat)
    # Do a little sanity checking:
    seq_lengths_list = []
    for record in alignment:
       seq_lengths_list.append(len(record))
    seq_lengths = set(seq_lengths_list)
    if len(seq_lengths) != 1:
        sys.stderr.write("Your alignment lengths aren't equal. Check your alignment file.")
        sys.exit(1)
    index = range(1, list(seq_lengths)[0]+1)
    return alignment, list(seq_lengths), index

def shannon_entropy(list_input):
    import math
    unique_base = set(list_input)
    M   =  len(list_input)
    entropy_list = []
    # Number of residues in column
    for base in unique_base:
        n_i = list_input.count(base) # Number of residues of type i
        P_i = n_i/float(M) # n_i (Number of residues of type i) / M(Number of residues in column)
        entropy_i = P_i*(math.log(P_i,2))
        entropy_list.append(entropy_i)
    sh_entropy = -(sum(entropy_list))
    return sh_entropy

def shannon_entropy_list_msa(alignment):
    # Calculate Shannon Entropy across the whole MSA
    shannon_entropy_list = []
    for col_no in range(len(list(alignment[0]))):
        list_input = list(alignment[:, col_no])
        shannon_entropy_list.append(shannon_entropy(list_input))
    return shannon_entropy_list

def main():
    alignment, seq_lengths, index = parseMSA(msa, alnformat, verbose)
    sel = shannon_entropy_list_msa(alignment)
    entropy_record = []
    for c1, c2 in zip(index, sel):
        record_to_append = c2
        entropy_record.append(record_to_append)
    return entropy_record

entropy_record = main()

# load entropy data and remove those at gap% beyond threshold

entropylist_S = pd.Series(entropy_record)

entropylist_S_gap_considered = entropylist_S.drop(gap_percent_bool_index_remove)
max_entropy_before_gaps = pd.Series.max(entropylist_S)
print('Maximum Shannon entropy alignment score before gap % considered: ', f'{max_entropy_before_gaps: .2f}')
max_entropy_after_gaps = pd.Series.max(entropylist_S_gap_considered)
print('Maximum Shannon entropy alignment score after gap % considered: ', f'{max_entropy_after_gaps: .2f}')

entropy_record_numpy = entropylist_S_gap_considered.to_numpy()
entropy_record_numpy.shape = (-1,len(entropylist_S_gap_considered))

print('Preparing sequences for comparison.')

# load sequences from alignment into list and control case

record_x_toward_seq_dataframe = []
sequence_records = []

for record_x in SeqIO.parse(input_sequence,"fasta"):
    record_x_toward_seq_dataframe = list(record_x.seq)
    record_x_toward_seq_dataframe_lower = [x.lower() for x in record_x_toward_seq_dataframe] 
    record_x_toward_seq_dataframe_ASCII = [ord(x) for x in record_x_toward_seq_dataframe_lower]
    sequence_records.append(record_x_toward_seq_dataframe_ASCII)

#generate dataframe of alignment from list

print('Generating sequence dataframe.')
sequence_dataframe = pd.DataFrame(sequence_records)

sequence_dataframe = sequence_dataframe.astype('int8')

# remove gapped positions

print('Removing gapped positions from analysis set.')
sequence_dataframe_gap_considered = sequence_dataframe.drop(gap_percent_bool_index_remove,axis=1)
print("Elapsed time: ~ " + str(int(time.time() - start_time)) + " seconds.")

# clear out unused items from memory

del sequence_records
del sequence_dataframe
gc.collect()

# prepare dataframe for storage of trial results (these columns are stripped away later if only one trial is performed)

record_sequence_trial_results['Total_trials'] = 0
record_sequence_trial_results['Outlier_instances'] = 0

# set trial counter
print('Beginning sequence trials.')
trial_count = 0

# avoid empty source dataframe

if depth_of_alignment//number_in_small_test == depth_of_alignment/number_in_small_test:
    times_to_sample_max_keep = (depth_of_alignment//number_in_small_test)
else:
    times_to_sample_max_keep = (depth_of_alignment//number_in_small_test) + 1

# define the calculation engine

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
                    # comparison_bool_series_match_index = np.where(comparison_bool_series_match==True) // in development for numba compilation
                    # comparison_bool_series_NOT_match_index = np.where(comparison_bool_series_NOT_match==True) // in development for numba compilation
                    # entropy_record_match = np.take(entropy_record_numpy,comparison_bool_series_match_index) // in development for numba compilation
                    # entropy_record_NOT_match = np.take(entropy_record_numpy,comparison_bool_series_NOT_match_index) // in development for numba compilation
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
    
    sequence_dataframe_gap_considered = sequence_dataframe_gap_considered.sample(frac=1) # shuffle master sequence dataframe
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
    
    # initiate numpy array for entropy calculation values
    
        entropy_array = np.empty((number_to_choose,number_to_choose),dtype=float)
        entropy_array[:] = np.nan
                          
    # calculations of match or not, and sum entropy values

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

# calculate interquartile range and outlier cutoff

        entropy_DF_analysis_values_list = entropy_DF_analysis.values.tolist()
        q25, q75 = np.nanpercentile(entropy_DF_analysis_values_list, 25), np.nanpercentile(entropy_DF_analysis_values_list, 75)
        iqr = q75 - q25

        CIQR = iqr * multiplier_on_interquartile_range
        lower_cutoff, upper_cutoff = q25 - CIQR, q75 + CIQR

# identify the outlier sequences using the interquartile range cutoff

        entropy_DF_analysis_above_cutoff = entropy_DF_analysis > upper_cutoff

        entropy_median_too_high = entropy_DF_analysis_above_cutoff.loc[entropy_DF_analysis_above_cutoff['Median'] == True]

        record_sequence_trial_results.loc[entropy_median_too_high.index,'Outlier_instances'] += 1
    print("Elapsed time: ~ " + str(int(time.time() - start_time)) + " seconds.")
    print("Estimated total time for analysis: ~ " + str(int(((time.time() - start_time))/(1+trial)*min_trials_for_each_sequence)) + " seconds.")

# print full distance matrix for analysis and generate a plot only if a single test of all sequences versus all sequences was performed

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
    
    # print figure describing trials
    
    #fig, ax = plt.subplots()
    #plt.figure(figsize=(8,8))
    #plt.xlabel('Sequence number within input file', fontsize=12)
    #plt.ylabel('Sequence number within input file', fontsize=12)
    #if depth_of_alignment < 2000:
    #    plotlabels = [i+1 for i in range(0,(depth_of_alignment-1),50)]
    #if depth_of_alignment >= 2000:
    #    plotlabels = [i+1 for i in range(0,(depth_of_alignment-1),500)]
    #plt.xticks(ticks=plotlabels,rotation=60)
    #plt.yticks(ticks=plotlabels)
    #plt.imshow(entropy_DF,cmap='cividis',interpolation='nearest')
    #plt.show
    #plt.colorbar()
    #plt.savefig(output_figure + '.pdf', dpi=600)

# prepare dataframe to generate FASTA files

record_sequence_trial_results['Fraction_positive'] = record_sequence_trial_results['Outlier_instances'] / record_sequence_trial_results['Total_trials']

record_seq_convert_to_string = []
for record in SeqIO.parse(input_sequence,"fasta"):
    record_seq_convert_to_string.append(str(record.seq))
    
acc_records_S = pd.Series(alignment_record_name_list)
sequence_records_S = pd.Series(record_seq_convert_to_string)
    
frame = { 'Accession': acc_records_S, 'Sequence': sequence_records_S }
FASTA_output_unclean_DF = pd.DataFrame(frame) 

# generating clean dataframes

if stringency_flag == 1: # minimal stringency
    FASTA_output_clean_DF = FASTA_output_unclean_DF.loc[record_sequence_trial_results['Fraction_positive'] != 1]
if stringency_flag == 2: # moderate stringency
    FASTA_output_clean_DF = FASTA_output_unclean_DF.loc[record_sequence_trial_results['Fraction_positive'] <= 0.5]
if stringency_flag == 3: # maximal stringency
    FASTA_output_clean_DF = FASTA_output_unclean_DF.loc[record_sequence_trial_results['Fraction_positive'] == 0]

# generating rejection dataframes

if stringency_flag == 1: # minimal stringency
    FASTA_output_reject_DF = FASTA_output_unclean_DF.loc[record_sequence_trial_results['Fraction_positive'] == 1]
if stringency_flag == 2: # moderate stringency
    FASTA_output_reject_DF = FASTA_output_unclean_DF.loc[record_sequence_trial_results['Fraction_positive'] > 0.5]
if stringency_flag == 3: # maximal stringency
    FASTA_output_reject_DF = FASTA_output_unclean_DF.loc[record_sequence_trial_results['Fraction_positive'] != 0]

# save clean FASTA file

print('Writing cleaned alignment as FASTA.')
print(FASTA_output_clean_DF)
FASTA_output_final_acc_list = FASTA_output_clean_DF.loc[:,'Accession'].values.tolist()
FASTA_output_final_seq_list = FASTA_output_clean_DF.loc[:,'Sequence'].values.tolist()

ofile = open(output_sequence, "w")
    
for seqi in range(len(FASTA_output_final_acc_list)):
    ofile.write(">" + FASTA_output_final_acc_list[seqi] + "\n" + FASTA_output_final_seq_list[seqi] + "\n")
ofile.close()

# save rejected FASTA file

print('Writing rejected sequences to FASTA.')
print(FASTA_output_reject_DF)
FASTA_output_rejected_acc_list = FASTA_output_reject_DF.loc[:,'Accession'].values.tolist()
FASTA_output_rejected_seq_list = FASTA_output_reject_DF.loc[:,'Sequence'].values.tolist()
    
ofile = open(output_rejected, "w")
for seqi in range(len(FASTA_output_rejected_acc_list)):
    ofile.write(">" + FASTA_output_rejected_acc_list[seqi] + "\n" + FASTA_output_rejected_seq_list[seqi] + "\n")
ofile.close()

# save analysis file as CSV

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

# provide runtime

print("Analysis complete.")
print("Total time for analysis: ~ " + str(int(((time.time() - start_time))/(1+trial)*min_trials_for_each_sequence)) + " seconds.")

# %%