#!/usr/bin/env python
import os
import sys
import shannon_support_read_finders as srf

if len(sys.argv) <5:
	print('need recon file')
	print('need shannon output dir')
	print('read.fasta output dir')
	print('is_check seq, 1 for check')
	print('')
	print('NOTE: output is named with comp_num')
	sys.exit()

recon_f = sys.argv[1]
shannon_out_dir = os.path.abspath(sys.argv[2])
read_fasta_dir = sys.argv[3]
is_check_seq = int(sys.argv[4]) == 1


if not os.path.exists(read_fasta_dir):
	os.makedirs(read_fasta_dir)

def get_comp_num_from_header(line):
	tokens = line.split()
	header = tokens[0]
	header_tokens = header.split('_')
	comp_num = int(header_tokens[1])
	return comp_num

def get_num_comp(dump_comp_dir):
	return len(os.listdir(dump_comp_dir))

dump_comp_dir = shannon_out_dir + '/' + 'components_dump_reads'
num_comp = get_num_comp(dump_comp_dir)

all_read_ids = set()

comp_fasta = {}

with open(recon_f) as f:
	for line in f:
		if line[0] == '>':
			header = line[1:]
			seq = next(f)
			if header[:6] == 'Contig':
				continue
			comp_num = get_comp_num_from_header(header)	
			if comp_num not in comp_fasta:
				comp_fasta[comp_num] = [(header, seq)]
			else:
				comp_fasta[comp_num].append((header, seq))


	
for comp_num, fasta in comp_fasta.items():
	dumped_read_file = shannon_out_dir + '/' + 'components_dump_reads' + '/' + 'comp' + str(comp_num)
	comp_read_list, comp_count_list, comp_prob_list = srf.get_all_reads(dumped_read_file)
	all_read_ids = set()
	
	comp_dir = read_fasta_dir + '/' + 'comp' + str(comp_num)
	if not os.path.exists(comp_dir):
		os.makedirs(comp_dir)
	
	for header_seq in fasta:
		header, seq = header_seq[0], header_seq[1]
		sh_header =  header.split()[0]

		read_ids_set, sum_rc = srf.fetch_read_ids_for_read(header, seq, shannon_out_dir, is_check_seq, comp_read_list, comp_count_list)
		comp_graph_out_fasta = comp_dir + '/' + str(sh_header) + '.read.fasta'
		srf.dump_comp_graph_reads(comp_graph_out_fasta, read_ids_set, comp_read_list, comp_count_list, comp_prob_list)

		all_read_ids.update(read_ids_set)

	all_reads_fasta = comp_dir + '/' + 'comp' + str(comp_num) + '.read.fasta'
	srf.dump_comp_graph_reads(all_reads_fasta, all_read_ids, comp_read_list, comp_count_list, comp_prob_list)


