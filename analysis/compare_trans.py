import sys
import tester
import os

os.chdir('./analysis')

python_path = 'python'

def run_cmd(s1):
	os.system(s1)

def compare_trans(reference,reconstr,reconstr_log):
	reconstr_per = reconstr_log + '_per'
	run_cmd(python_path + ' ' +  'parallel_blat_python.py ' + reconstr + ' ' + reference + ' ' + reconstr_per)
	tester.analyzer_blat_noExp(reconstr_per,reconstr_log,'  ',1000000)
	#if 0: #false_positive:
	#	tester.false_positive(reconstr,reconstr_per,reconstr_rev_log)

reference = sys.argv[1]
reconstr = sys.argv[2]
reconstr_log = sys.argv[3]
#print(reference)
#print(reconstr)
#print(reconstr_log)
compare_trans(reference,reconstr,reconstr_log)
