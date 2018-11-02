#!/usr/bin/env	python
import os
import sys
import subprocess
import shannon_support_read_finders as srf

if len(sys.argv) != 4:
	print('need input a fasta file, currently only support s_***, comp_***')
	print('need shannon output dir')
	print('enter non-zero to check if reads actually match exactly to the seq, default to 0(false)')
	sys.exit()

input_fasta = sys.argv[1]
shannon_out_dir = os.path.abspath(sys.argv[2])
is_check_seq = sys.argv[3] != '0'

with open(input_fasta) as f:
	for line in f:
		if line[0] == '>':
			header_line = line[1:]
			sh_header = header_line.split()[0]
			comp_num = srf.get_comp_num(header_line)
			#sum_rc = srf.get_transcript_rc(header_line)
			
			seq = next(f).split()[0]

			dumped_reads_path = shannon_out_dir + '/' + 'components_dump_reads' + '/' + 'comp' + str(comp_num)
			comp_read_list, comp_count_list, comp_prob_list = srf.get_all_reads(dumped_reads_path)
			read_ids_set, sum_rc = srf.fetch_read_ids_for_read(header_line, seq, shannon_out_dir, is_check_seq, comp_read_list, comp_count_list)
			print('num_report_from_sh ', sum_rc)
			
			#writer.write('>' + header + '\n')
			comp_graph_fasta = sh_header + '.read.fasta'
			num_supportive_read = srf.dump_comp_graph_reads(comp_graph_fasta, read_ids_set, comp_read_list, comp_count_list, comp_prob_list)
			print('num_supportive_read', num_supportive_read)
			print('num uniq read id   ', len(read_ids_set))

sys.exit(0)
