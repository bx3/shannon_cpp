import sys

node_prefix = 'Graph_nodes_'
path_prefix = 'Graph_paths_'
read2_prefix = 'Graph_reads2_'

def check_seq_align(seq, read_seq):
	if seq.find(read_seq) == -1:
		print('find a non-match read')
		print('read: ' + read_seq)
		print('fasta:' + seq)
		return False
			
	return True

def get_comp_num(line):
	tokens = line.split()
	header = tokens[0]
	header_tokens = header.split('_')
	comp_num = int(header_tokens[1])
	return comp_num

def get_transcript_rc(line):
	tokens = line.split()
	return int(get_single_content_after_colon(tokens[2]))

def dump_comp_graph_reads(comp_graph_fasta, read_ids_set, comp_read_list, comp_count_list, comp_prob_list):
	num_supportive_read = 0
	with open(comp_graph_fasta, 'w') as writer:
		for read_id in read_ids_set:
			read_seq = comp_read_list[read_id]
			count = comp_count_list[read_id]
			prob = comp_prob_list[read_id]
			num_supportive_read += count
			#effective_rc = int(count/prob)
			writer.write('@' + str(count)+ '\n')#'\t'+str(prob)+'\t'+str(effective_rc) + '\n')
			writer.write(read_seq + '\n')

	return num_supportive_read
	#print('num_supportive_read', num_supportive_read)
	#print('num uniq read id   ', len(read_ids_set))

def get_content_after_colon(line, delimiter=','):
	colon_index = line.find(':')
	if colon_index == len(line)-1:
		return []
	content_list = line[colon_index+1:].split(delimiter)
	return content_list
def get_single_content_after_colon(line):
	colon_index = line.find(':')
	content = line[colon_index+1:].split()
	return content[0] 

def parse_header_info(line):
	tokens = line.split()
	header = tokens[0]
	
	num_sf = int(get_single_content_after_colon(tokens[1]))
	sum_rc = int(get_single_content_after_colon(tokens[2]))
	rc_1 = int(get_single_content_after_colon(tokens[3]))
	rc_2 = int(get_single_content_after_colon(tokens[4]))
	rc_3 = int(get_single_content_after_colon(tokens[5]))
	vd_path_raw = tokens[6]
	
	path2_ids_raw = tokens[7] # last oone is ,
	if path2_ids_raw[-1] == ',':
		path2_ids_raw = path2_ids_raw[:-1]
	path3_ids_raw = tokens[8]	# last oone is ,
	if path3_ids_raw[-1] == ',':
		path3_ids_raw = path3_ids_raw[:-1]

	vd_path = get_content_after_colon(vd_path_raw, '->')
	path2_ids = get_content_after_colon(path2_ids_raw)
	path3_ids = get_content_after_colon(path3_ids_raw)
	header_tokens = header.split('_')
	comp_num = int(header_tokens[1])
	graph_num = int(header_tokens[3])
	return header, comp_num, graph_num, vd_path, path2_ids, path3_ids, sum_rc
	
	
def parse_single_node_header_info(line):
	tokens = line.split()
	#print(tokens)
	header = tokens[0]
	header_tokens = header.split('_')
	comp_num = int(header_tokens[1])
	
	vd_rc = get_single_content_after_colon(tokens[1])
	i = vd_rc.find('(')
	j = vd_rc.find(')')
	rc = int(vd_rc[i+1:j])
	return header, comp_num, vd_rc, rc


#def get_Path2_supportive_reads_id():
def get_files(comp_num, graph_num, shannon_out_dir):
	node_path = shannon_out_dir + '/' + 'comp_graph' + '/' + node_prefix + str(comp_num) + '/' + 'node' + str(graph_num)
	
	path_path = shannon_out_dir + '/' + 'comp_graph' + '/' + path_prefix + str(comp_num) + '/' + 'path' + str(graph_num)

	read_path = shannon_out_dir + '/' + 'comp_graph' + '/' + read2_prefix + str(comp_num) + '/' + 'read' + str(graph_num)
	return node_path, path_path, read_path

def get_node_id_and_read_count(vd):
	i = vd.find('(')
	j = vd.find(')')
	vd_id  = vd[:i]
	support_read_count = vd[i+1:j]
	return vd_id, support_read_count
	
def get_read_id_and_read_count(read_str):
	i = read_str.find('(')
	j = read_str.find(')')
	read_id  = int(read_str[:i])
	support_read_count = int(read_str[i+1:j])
	return read_id, support_read_count


def get_read_ids_for_node_file(header, node_file, vd_path, read_list, count_list, fasta_seq, is_check_seq):
	vd_id_rc = {}
	for vd in vd_path:
		#print('vd', vd) 
		vd_id, rc = get_node_id_and_read_count(vd)
		vd_id_rc[vd_id] = rc
	supportive_read_ids = set()
	#print('node_file', node_file)
	#print('vd_id_rc', vd_id_rc)
	with open(node_file) as f:
		next(f)	
		for line in f:
			tokens = line.split()
			vd_id = tokens[0]
			if vd_id in vd_id_rc:
				#print(tokens)
				bases = tokens[1]
				count = tokens[2]
				supportive_read_count = tokens[3]
				reads_ID = tokens[4:]
				for read_str in reads_ID:
					read_id, rc = get_read_id_and_read_count(read_str)
					if is_check_seq:
						if not check_seq_align(fasta_seq, read_list[read_id]):
							print('header', header)
							print('from node, read id not align', read_id)
							print('node id', vd_id)
							sys.exit()
							
						supportive_read_ids.add(read_id)
	
	return list(supportive_read_ids)

def get_read_ids_for_path_file(header, path_file, path3_ids_str, read_list, count_list, fasta_seq, is_check_seq):
	if len(path3_ids_str) == 0:
		return []
	i = 0;
	j = 0
	read_ids = []
	path3_ids = [int(k) for k in path3_ids_str]
	path_ids_set = set(path3_ids)
	max_line = max(path3_ids)
	supportive_read_ids =set() 

	with open(path_file) as f:
		next(f)
		for line in f:
			if j >  max_line:
				break
			if j in path_ids_set:		
				tokens = line.split()
				num_node = int(tokens[0])
				vd_list = tokens[1:1+num_node]
				read_ids = tokens[1+num_node+1:]
				for read_str in read_ids:
					read_id, rc = get_read_id_and_read_count(read_str)
					if is_check_seq and not check_seq_align(fasta_seq, read_list[read_id]):
						print('header', header)
						print('from path, read id not align', read_id)
						print('path id', j)
						sys.exit()
					supportive_read_ids.add(read_id)
				
				i += 1
			j += 1
	return list(supportive_read_ids)

def get_read_ids_for_read_file(header, read_file, path2_ids_str, read_list, count_list, fasta_seq, is_check_seq):
	if len(path2_ids_str) == 0:
		return []
	i = 0;
	j = 0
	read_ids = []
	path2_ids = [int(k) for k in path2_ids_str]
	path_ids_set = set(path2_ids)
	max_line = max(path2_ids)
	supportive_read_ids =set() 

	with open(read_file) as f:
		next(f)
		for line in f:
			if j > max_line:
				break
			if j in path_ids_set:		
				tokens = line.split()
				vd1 = tokens[0]
				vd2 = tokens[1]
				read_count = int(tokens[2])

				read_ids = tokens[3:]
				for read_str in read_ids:
					read_id, rc = get_read_id_and_read_count(read_str)
					if is_check_seq and not check_seq_align(fasta_seq, read_list[read_id]):
						print('header', header)
						print('from read, read id not align', read_id)
						print('read line', j)
						sys.exit()

					supportive_read_ids.add(read_id)

			j += 1
	return list(supportive_read_ids)

def get_all_reads(dumped_read_file):
	reads_list  = []
	count_list = []
	prob_list = []
	with open(dumped_read_file) as f:
		for line in f:
			if line[0] == '>':
				tokens = line.split()
				read_id = int(tokens[0][4:])
				read_count = int(tokens[1])
				prob = float(tokens[2])
				bases = next(f)
				bases = bases.split()[0]
			
				reads_list.append(bases) 
				count_list.append(read_count)
				prob_list.append(prob)
			else:
				print('dump fasta not align to node')
				print('fasta ' + bases)
				print('node  ' + seq)
				sys.exit()
	
	return reads_list, count_list, prob_list

	
def fetch_read_ids_for_read(line, seq, shannon_out_dir, is_check_seq, read_list, count_list):
	comp_num = -1
	all_read_ids = []
	header = ''
	sum_rc = 0
	
	if line[0] == 's':
		is_single_node = True
		header, comp_num, vd_rc, sum_rc = parse_single_node_header_info(line)		

		node_path = shannon_out_dir + '/' + 'comp_graph' + '/' + node_prefix + str(comp_num) + '/' + 'node_s'
		all_read_ids = get_read_ids_for_node_file(header, node_path, [vd_rc], read_list, count_list, seq, is_check_seq)

	elif line[:4] == 'comp':
		header, comp_num, graph_num, vd_path, path2_ids, path3_ids, sum_rc = parse_header_info(line)

		node_path, path_path, read_path = get_files(comp_num, graph_num, shannon_out_dir)
		nodes_read_ids = list(set(get_read_ids_for_node_file(header, node_path, vd_path, read_list, count_list, seq, is_check_seq)))
		reads_read_ids = list(set(get_read_ids_for_read_file(header, read_path, path2_ids, read_list, count_list, seq, is_check_seq)))
		paths_read_ids = list(set(get_read_ids_for_path_file(header, path_path, path3_ids, read_list, count_list, seq, is_check_seq)))
		all_read_ids = set(nodes_read_ids + reads_read_ids + paths_read_ids)
	return all_read_ids, sum_rc
