#!/usr/bin/env python
import os
import sys
import subprocess

# setting
single_read_path = '../sample/SE_read.fasta'
paired_read_1_path = '../sample/PE_read_1.fasta'
paired_read_2_path = '../sample/PE_read_2.fasta'

 
output_dir = '../sample/sample_output'
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

single_output_dir = '../sample/sample_output/single_out'
paired_output_dir = '../sample/sample_output/paired_out'

single_reconstruct_fasta_path = single_output_dir + '/' + 'reconstructed_seq.fasta'
paired_reconstruct_fasta_path = paired_output_dir + '/' + 'reconstructed_seq.fasta'

single_length = 50
paired_1_length = 100
paired_2_length = 100

num_threads = '4'
sampling = '0'
max_memory = '4G'
sort_threads = '1'
random_seed = '31'
min_output_fasta_len = '80'

def run_single_reads_Shannon():
	cmd = ['../Shannon_RNASeq_Cpp', 'shannon',
								   '-l', str(single_length),
								   '-s', single_read_path,
								   '-o', single_output_dir,
								   '-t', num_threads,
								   '-g', sampling,
								   '-m', max_memory,
								   '-u', sort_threads,
								   '-e', min_output_fasta_len,
								   '--random_seed', random_seed]
	#print('run single end reads Shannon command')
	#print(' '.join(x for x in cmd))
	result = os.system(' '.join(x for x in cmd))
	if result != 0:
		print('Shannon single reads input returns error code', result)
		sys.exit()

	return 0


def run_paired_reads_Shannon():
	length_arg = str(paired_1_length) + ' ' + str(paired_2_length)
	print(length_arg)
	cmd = ['../Shannon_RNASeq_Cpp', 'shannon', 
					  '-i', str(paired_1_length) + ' ' + str(paired_2_length),
					  '-p', paired_read_1_path + ' ' + paired_read_2_path,
					  '-o', paired_output_dir,
					  '-t', num_threads,
					  '-g', sampling,
					  '-m', max_memory,
					  '-u', sort_threads,
					  '-e', min_output_fasta_len,
					  '--random_seed', random_seed]
	#print('run paired end reads Shannon command')
	#print(' '.join(x for x in cmd))
	result = os.system(' '.join(x for x in cmd))
	if result != 0:
		print('Shannon paired reads input returns error code', result)
		sys.exit()

	return 0

def check_if_file_nonempty(f):
	if(not os.path.exists(f)):
		print('ERROR: Shannon Installation failed, it failed to create output ', f)
		sys.exit()
	else:
		if(os.stat(f).st_size == 0):
			print('ERROR: Shannon Installation failed, it created empty output ', f)
			sys.exit()
		else:
			print('Shannon reconstructed transcriptomes for', f)


#def check_if_reproduce_same_output(new_f, seed_31_f):



run_single_reads_Shannon()

run_paired_reads_Shannon()

check_if_file_nonempty(single_reconstruct_fasta_path)
print('************************************')
print('Shannon passes: for single end reads')
print('************************************')

check_if_file_nonempty(paired_reconstruct_fasta_path)
print('************************************')
print('Shannon passes: for paired end reads')
print('************************************')
