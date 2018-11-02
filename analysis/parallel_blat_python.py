import subprocess
import os
import pdb
import math
import run_parallel_cmds

def cut_file(in_name,out_name,line_start,line_end):
	os.system('gawk \'NR > ' + str(line_end) + ' { exit } NR >= ' + str(line_start) +  '\' '+ in_name + ' > ' + out_name )

def parallel_blat(target_fasta,query_fasta,out_file,QUERY_SPLIT,nJobs=60):
	'''Function takes in target,query and output file. parallelizes blat by running GNU parallel
	- Currently only parallelizes on query space
	- Also assumes that query fasta file takes two lines per sequence (not wrapped)'''
	target_length = float(subprocess.check_output('grep -c \'>\' ' + target_fasta,shell=True))
	query_length = float(subprocess.check_output('grep -c \'>\' ' + query_fasta,shell=True))
	os.system( 'awk \'/^>/{print s? s"\\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}\' ' +query_fasta + ' > '+ query_fasta +'_nospace')
	#os.system( 'awk \'/^>/{print s? s"\\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}\' ' +target_fasta + ' > '+ target_fasta)
	query_fasta = query_fasta + '_nospace'
	#TARGET_SPLIT = 1
	#QUERY_SPLIT = 4
	#Alernately
	#QUERY_SPLIT = min(int(math.ceil(float(query_length)/float(target_length))),50)
	#QUERY_SPLIT = max(int(math.ceil(float(query_length)/float(target_length))),500) 
	#QUERY_SPLIT = int(min(QUERY_SPLIT,query_length))
	#QUERY_SPLIT= 100
	#pdb.set_trace()
	print('Query length: ' +str(query_length) + ' Target length: ' + str(target_length) + ' Query Split: ' + str(QUERY_SPLIT))
	split_size = int(math.floor(float(query_length)/QUERY_SPLIT))
	'''if split_size % 2 !=0:
		split_size +=1'''
	'''if query_length <= float(target_length):
		print('Cannot parallelize on query. Running Vanilla Blat')
		os.system('./blat -noHead ' + target_fasta + ' ' +  query_fasta + ' ' + out_file)	
		return'''
	query_size = int(query_length / QUERY_SPLIT)
	if query_size%2 != 0:
		query_size += 1

	split_dir = out_file[:-1] + '_split'
	if not os.path.exists(split_dir):
		os.makedirs(split_dir)
	subprocess.call(['split', '-dl', str(query_size), query_fasta, split_dir + '/split'])
	n = 0
	cmds = []
	split_files = os.listdir(split_dir)
	for f in split_files:
		cmds.append('./blat -noHead '+ target_fasta + ' ' + split_dir+'/'+f + ' ' + out_file + '_' + str(n))
		#print('blat -noHead '+ target_fasta + ' ' + split_dir+'/'+f + ' ' +out_file + '_' + str(n))
		n += 1

	cmds = tuple(cmds)
	run_parallel_cmds.run_cmds(cmds,nJobs)

	#print('parallel blat -noHead ' + target_fasta + ' ' + query_fasta + '_{} ' +out_file + '_{} ::: ' + q_str )
	#os.system('time parallel blat -noHead ' + target_fasta + ' ' + query_fasta + '_{} ' +out_file + '_{} ::: ' + q_str  )
	#os.system('parallel blat ' + target_fasta + ' ' + query_fasta + '_{} ' +out_file + '_{} ::: {1..' + str(QUERY_SPLIT) + '}' )
	#os.system('sort -k 10 ' + out_file + '_* > ' + out_file)
	os.system('cat ' + out_file + '_* > ' + out_file)
	os.system('rm ' + split_dir + '/*')
	os.system('rm ' + out_file + '_*' )
	os.system('rm ' + query_fasta + '_*' )

def main():
    import sys
    args = sys.argv[1:]
    nJobs=60
    if '--nJobs' in args:
    	ind1 = args.index('--nJobs')
    	nJobs = int(args[ind1+1])
    	args = args[:ind1]+args[ind1+1:]    
    Query_split = 100

    
    parallel_blat(args[0],args[1],args[2],Query_split, nJobs)


if __name__ == '__main__':
    main()

