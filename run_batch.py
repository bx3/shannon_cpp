#!/usr/bin/env python
import sys
import os
import subprocess
import multiprocessing
import argparse

parser = argparse.ArgumentParser()
# Single
parser.add_argument("-s", dest='single_path')
parser.add_argument("-l", dest='single_length')

# Paired
parser.add_argument("-1", dest='left_path')
parser.add_argument("-2", dest='right_path')
parser.add_argument("-i", nargs='*', dest='paired_lengths')
parser.add_argument("-o", dest='output_dir', required=True)

def run_process(cmd):
    print(cmd)
    os.system(cmd)

def run_cmds(cmds,noJobs):
    #cmds is a tuple of strings
    #noJobs is a an integer with the no of Jobs to run in parallel
    p = multiprocessing.Pool(noJobs)
    p.map(run_process,cmds)

def get_inputs_cmd(args):
    args = ' '.join(args)
    inputs_cmd = []
    cmd = ''
    for c in args:
        if ',' == c:
            if cmd != '':
                t = cmd.strip().split(' ')
                inputs_cmd.append(t)
                cmd = ''
        else:
            cmd+=str(c)

    if len(cmd) != 0:
        inputs_cmd.append(cmd.strip().split(' ') )
    return inputs_cmd


def delete_arg(args, label):
    index = args.index(label)
    if '-' in args[index+1]:
        print("detect args format error -x followed by -x. Your input ", args)
    del args[index: index+2]

def find_and_delete(args, label, shannon_cmd):
    count = args.count(label)

    if count == 0:
        return -1
    else:
        index = args.index(label)
        if '-' in args[index+1]:
            print("detect args format error -x followed by -x. Your input ", args)
        value = args[index+1]
        if ',' in value:
            value = value.replace(',', '')

        shannon_cmd.append(label)
        shannon_cmd.append(value)
        del args[index: index+2]
        return value
        
def find_value(args, label):
    count = args.count(label)
    if count == 0:
        return -1
    else:
        index = args.index(label)
        if '-' in args[index+1]:
            print("detect args format error -x followed by -x. Your input ", args)
        value = args[index+1]
        if ',' in value:
            value = value.replace(',', '')
        del args[index: index+2]
        return value

def parse_global(args, global_cmd):
    jobs =  find_value(args, '-n')
    threads = find_and_delete(args, '-t', global_cmd)
    kmer_length =  find_and_delete(args, '-k', global_cmd)
    sampling = find_and_delete(args, '-g', global_cmd)
    memory = find_and_delete(args, '-m', global_cmd)
    sort_thread = find_and_delete(args, '-u', global_cmd)
    min_length = find_and_delete(args, '-e', global_cmd)

    if jobs == -1 or threads==-1:
        print('need ')
        print('jobs: -n')
        print('threads: -t')
        sys.exit(0)
    
    return jobs, threads, kmer_length, sampling, memory, sort_thread, min_length, args 

def print_cmd_help():
    print("")
    print("This is a batch processing helper for handling multiple input seperated by comma")
    print("There are Two types of argument options")
    print("")
    print("LOCAL arguments relates to data, used for specifying input, output path and length")
    print("LOCAL arguments are seperated by comma")
    print("-s   single ended read path")
    print("-l   single ended max read length")
    print("-1   left input path")
    print("-2   right input path")
    print("-i   two values, the first for left input max length, the second for right input max length")
    print("-o   output directory")
    print("")
    print("Global arguments controls software processing scheme, applys to all inputs")
    print("-n   number of processes running concurrently, where each process handles one input")
    print("-k   kmer length")
    print("-t   number of threads per process")
    print("-g   given that many reads are redundant,")
    print("     setting g to subsample matching reads")
    print("     to speedup later works, where g can be")
    print("     interpreted as coverage depth. g=0: no")
    print("     subsample; g>0: probabilistic sample;")
    print("     g<0: sample every (-g) read")
    print("-u   number of threads used by linux sort program")
    print("-m   specify amount of memory to be used in multibridge step, ignoring -m leads to using all memory")
    print("-e   minimal reconstructed output length")
    print("")
    print("EXAMPLE of two inputs. The first inputs(A) contain both single and paired pair, the second inputs(B) contain only single read")
    print("./run_batch.py cmd -n 2 -t 12 ,-s path_to_single_read-A -l max_length_of_single_read-A -1 path_to_left_read-A -2 path_to_right_read-A  -i left_read_max_length-A right_read_max_length-A -o outdir-A, -s path_to_single_read-B -l max_length_of_single_read-B -o outdir-B")
    print("")

def print_file_help():
    print("")
    print("-n   number of processes running concurrently, where each process handles one input")
    print("-f   input file, where each line of the file is a valid shannon_cpp command, this python helper file simply run each line of commands included in the file ")
    print("")


def get_shannon_cmd(tokens):
    shannon_cmd = ['shannon_cpp', 'shannon']
    tokens_str = ' '.join(tokens)
    args = parser.parse_args(tokens)

    single_path = args.single_path
    single_length = args.single_length
    left_path = args.left_path
    right_path = args.right_path
    paired_lengths = args.paired_lengths
    output_dir = args.output_dir

    if single_path == None and (left_path==None or right_path==None):
        print("require input files. Your input is ", tokens_str)
        sys.exit(0)

    if single_path != None:
        if single_length == None:
            print("need input max length. Your input is ", tokens_str)
            sys.exit(0)

        shannon_cmd.append('-s')
        shannon_cmd.append(single_path)
        shannon_cmd.append('-l')
        shannon_cmd.append(single_length)

    if left_path != None and right_path != None:
        if paired_lengths == None or len(paired_lengths) != 2:
            print("needs a pair of length. Your input is ", tokens_str)
            sys.exit(0)
        left_length = paired_lengths[0]
        right_length = paired_lengths[1]

        shannon_cmd.append('-p')
        shannon_cmd.append(left_path)
        shannon_cmd.append(right_path)
        shannon_cmd.append('-i')
        shannon_cmd.append(left_length)
        shannon_cmd.append(right_length)

    elif (left_path == None and right_path != None) or (left_path != None and right_path == None):
        print("need to supply both left and right path. Your input is ", tokens_str)
        sys.exit(0)

    shannon_cmd.append('-o')
    shannon_cmd.append(output_dir)

    return shannon_cmd



if len(sys.argv) < 2 or (sys.argv[1]!='file' and sys.argv[1]!='cmd'):
    print("need to choose an input method:")
    print("command line: cmd")
    print("input file: file")
    print("use --help after choosing input methods")
    sys.exit(0)

choice = sys.argv[1]
args = sys.argv[2:]

shannon_cmds = []
jobs = -1

if choice == 'cmd':
    if args.count('--help') > 0 or args.count('-h') > 0:
        print_cmd_help()
        sys.exit(0)

    global_cmd_list = []
    jobs, threads, kmer_length, sampling, memory, sort_thread, min_length, remain_args = parse_global(args, global_cmd_list)
    inputs_cmd = get_inputs_cmd(remain_args)
    if len(inputs_cmd) == 0:
        print("need to supply input -s -l ... Your input ", args)
        parser.parse_args('')
        sys.exit(0)

    for in_cmd in inputs_cmd:
        shannon_cmd = get_shannon_cmd(in_cmd) 
        shannon_cmd += global_cmd_list
        shannon_cmd_str = ' '.join(shannon_cmd)
        shannon_cmds.append(shannon_cmd_str)

if choice == 'file':
    if args.count('--help') > 0 or args.count('-h') > 0:
        print_file_help()
        sys.exit(0)

    file_cmd = []
    jobs = find_value(args, '-j', file_cmd)
    filename = find_and_delete(args, '-f', file_cmd)
    if jobs == -1:
        print("need to specify number of process")
        sys.exit(0)
    if filename == -1:
        print("need to specify input filename")
        sys.exit(0)

    with open(filename) as f:
        for line in f:
            tokens = line.split()
            shannon_cmd_str = ' '.join(tokens)
            shannon_cmds.append(shannon_cmd_str)

jobs = int(jobs)
cmds = tuple(shannon_cmds)
print(cmds)
run_cmds(cmds, jobs)

