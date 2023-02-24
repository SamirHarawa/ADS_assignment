#!/usr/bin/env python3

#***************************************************************************************************************************************#
#                                                                                                                                       #
#       Script name:    task_1.py                                                                                                       #
#       Description:    Counts k-mers for a given string and length of k-mer.                                                           #
#                       Uses suffix array.                                                                                              #
#                       Uses  "timeit" module to capture algorithm execution.                                                           #                     
#       Module:         Algorithms and Data Structures in Bioinformatics - ADS 612                                                      #
#       Submitted By:   Vita Nyasulu, Samir Adrian & Limbani Thengo                                                                     #
#       Date:           24 February 2023                                                                                                #
#                                                                                                                                       #
#***************************************************************************************************************************************#

from os import path

"""
Accepts a FASTA file and returns a concatenated string of the passed DNA sequence
"""
def get_sequence_from_fasta(path_to_fasta):
    # Read sequence file and handle common file opening Exceptions
    try:
        # Using with so that we don't worry about closing the file after opening
        with open(path_to_fasta) as sequence_file:
            sequence = sequence_file.read()

    except FileNotFoundError:
        print(path_to_fasta + r" not found")   

        # Stop execution
        print("Stopping script")
        exit() 

    except NameError:
        print(path_to_fasta + r" not found") 

        # Stop execution
        print(r"Stopping script")
        exit()

    else:
        # Get number of lines so that we don't work with an empty FASTA file
        num_lines = len(sequence.split('\n'))

        # one line sequence 
        concatenated_sequence = ""

        # More than one line concatenate
        if(num_lines > 0):
            for line in sequence.splitlines():
                if line.startswith(r'>'):
                    # Header line skip
                    pass
                else:
                    concatenated_sequence += line.strip()
        else:
            concatenated_sequence += sequence
    
    # Return concatenated sequence
    return concatenated_sequence
    
"""
Accepts a sequence string / text
Returns a suffix array list for passed string
"""  
def get_suffix_array(sequence):
    text = sequence
    
    # Suffix array for sorted suffixes
    suffix_array = []

    for index in range(len(sequence)): 
        suffix_array.append([index, text[index:]])
    
    # Sort "suffix array" list by suffix value
    suffix_array.sort(key = lambda x : x[1])

    return suffix_array



"""
Accepts string an array of suffixes for a string and a second parameter for k-mer
Prints some information on running time, length of passed string
"""
def count_kmers(suffix_array, k):
    # kmer list and their counts
    kmer_counts = {}

    # Length of suffix array
    len_suffix_arr = len(suffix_array)

    kmer = k

    # Current suffix 
    curr_prefix = r""

    # Previous suffix
    prev_prefix = r""

    # Loop through suffix array and compare suffixes from previous one to the next if they match 
    for i in range(len_suffix_arr):
        # Get current suffix 
        suffix = suffix_array[i][1]
        curr_prefix = suffix_array[i][1][0 : kmer]
        
        if i > 0:
            prev_prefix = suffix_array[i - 1][1][0 : kmer]
        else:
            prev_prefix = suffix_array[i][1][0 : kmer]

        # Only evaluate current suffix if it is >= kmer
        if len(suffix) >= k:
            # Check if we already have this prefix in our kmers dictionary
            if curr_prefix in kmer_counts:

                # Check if the prefix of size kmer for this suffix is equal to the prefix of size kmer of the previous suffix in the array list
                if curr_prefix == prev_prefix:
                    kmer_counts[curr_prefix] = kmer_counts[curr_prefix] + 1
            else: 
                # if not there already add it still with a count of 1
                kmer_counts[curr_prefix] = 1
            
    # Return dictionary of kmers (keys) and their counts
    return kmer_counts



# Execute statements below only as a script (will not execute if imported as a module)
if __name__ == "__main__":

    # Fasta file to work with
    sequence_file = r'wuhan.fasta'

    # data directory in user's HOME directory
    data_directory = r'nucleotides'

    # User's HOME directory
    user_home = path.expanduser('~')

    # Full path to sequence Fasta file 
    path_to_fasta = path.join(user_home, path.join(data_directory, sequence_file))

    # k-mer 
    kmer = 3

    # Sequence 
    sequence = get_sequence_from_fasta(path_to_fasta)

    kmers_dict = count_kmers(get_suffix_array(sequence), kmer)

    for key in kmers_dict:
        print(str(key) + r" occurs " + str(kmers_dict[key]) + r" times")