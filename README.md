# Shannon C++ RNA sequencing

This is a **C++** implementatin for the work **ShannonRNA: an Information-Optimal de Novo RNA-Seq Assembler**.  A **python** implementation can be found at *http://sreeramkannan.github.io/Shannon/*, which served as a foundation for the current **C++** implementation. However, this C++ version significantly improves both time and memory efficiency, roughly 10x, so you would prefer using this one.

## Content
- [Prerequisites](#prerequisites)
- [Installing](#installing)
- [Input requirement](#input-requirement)
- [Getting started](#getting-started)
- [Flowchart](#flowchart)
- [Output file structure](#output-file-structure)
- [Usage](#usage)
- [Config file](#config-file)
- [Memory parameter](#memory-parameter)
- [Versioning](#versioning)
- [Author](#author)
- [License](#license)
- [Acknowledgments](#acknowledgments)
### Prerequisites

* ubuntu (tested on 15.10, or higher)
* g++ (tested on 5.5.0, or higher, g++ 4.8.5 tested having issue to link boost.program_options)

### Installing
* install python 
    * sudo apt install python2.7 python-pip 
* Install METIS
	* sudo apt-get install libmetis-dev
	* OR follow instructions from [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/download) .
* Install cmake 
	* sudo apt-get install cmake
* Install Jellyfish
	* follow instrctions from [JellyFish](http://www.cbcb.umd.edu/software/jellyfish/) 
	* tested on 2.2.9, but any higher 2.x.x version would work
* Install boost library (May be the most painful part if you have an old one on system, cmake might not link properly)
	* follow instruction from [boost](http://www.boost.org/doc/libs/1_66_0/more/getting_started/unix-variants.html)
	* ~~sudo apt-get install libboost-all-dev~~, since it would install an older version, which would not compile
	* remember to install the non-header part of boost. This implementation needs
		* Boost.ProgramOptions
		* Boost.Filesystem 
		* Boost.System
* Install glpk 	
	*  install version 4.62 or higher, follow instructions from [GLPK](https://en.wikibooks.org/wiki/GLPK/Linux_OS)
	*  ~~sudo apt-get install libglpk-dev~~, since it would install an older version, which gives error
* Install google sparsehash 
	* follow instructions from [google sparsehash](https://github.com/sparsehash/sparsehash)
* hopscotch -- https://github.com/Tessil/hopscotch-map (already included in the repo, no action needed)
* Download sparsepp -- https://github.com/greg7mdp/sparsepp
	* download from link above, extract and place the whole files in the directory where CMakeLists.txt resides
* Syrupy -- https://github.com/jeetsukumaran/Syrupy
    * used for memory profiling
    * No action needed, already included

* After all, go to CMakelist.txt directory and type
```
    cmake .
    make
```
	
## Input requirement
* For noisy input reads, Rcorrector is recommanded prior to Shannon processing
    https://github.com/mourisl/Rcorrector
* Input can take a single-ended FASTA file or a pair-ended FASTA file or BOTH.
* all reads needs to be single-lined fasta file (multiline not supported)
```
Should look like this
>ABC
AGTGCGGCTTTCGTATTTGCTGCTCGTGTTTACTCTCACAAACTTGACCTGCACGCCAAAGAGATGCTTCTTGTGGAACTCGACA
>ABC2
AGTGCGGCTTTCGTATTTGCTGCTCGTGTTTACTCTCACAAACTTGACCTGCACGCCAAAGAGATGCTTCTTGTGGAACTCGACA
...
```
```
Instead of
>ABC
AGTGCGGCTTTCGTATTTGCTGCTCGTGTTTACTCTC
ACAAACTTGACCTGCACGCCAAAGAGATGCTTCTTGT
GGAACTCGACA
>ABC2
AGTGCGGCTTTCGTATTTGCTGCTCGTGTTTACTCTC
ACAAACTTGACCTGCACGCCAAAGAGATGCTTCTTGT
GGAACTCGACA
...
```
* Some useful command for converting fastq to fasta at https://gif.biotech.iastate.edu/converting-fastq-fasta


## Getting Started
* Download and install according to Installing section above. After installation, compile and run executable **Shannon_D_RNA_seq** with following options
```
Usage: Shannon_D_RNA_seq <CMD> arguments

shannon             run entire Shannon procedures
partition           run component partitions algorithm
multi-graph         run multibridge or sparse-flow or both on all components 
multi-bridge        run multibridge algorithm for one component
sparse-flow         run sparse-flow algorithm for one component
find-rep            Reduce redundant reconstructed transcripts after sparse-flow
ref-align           Align to reference file, and generate an evaulation file
custom              use json configuration file

Usage: Shannon_D_RNA_seq <CMD> --help to see options
```
### CMD\<shannon>
Shannon performs denovo transcriptome assembly by first partitioning reads into components, which can be individually run with CMD\<partition>. Shannon then builds a de-bruijn graph for each component and runs CMD\<multi-bridge> and CMD\<sparse-flow> on that graph. To process all components at the same time, CMD\<multi-graph> spawns multiple processes to facilitate parallel computing. At this step, Shannon has obtained reconstructed transcriptome and it uses CMD\<find-rep> to remove redundancy. After all, CMD\<ref-align> is optional, and allows comparison to the *ground true* result using BLAT. 

CMD\<shannon> takes raw reads (or Rcorrected reads) as its input.

**-s** specifies a single-ended fasta file at **path_to_single_ended_read** 
**-l** specifies the maximal length among all read, its main purpose is for accessing read efficiently since all reads are stored in array, -l has indexing role. Any reads with higher length are truncated, and the rest of the reads don't lose any information. If input reads have been preprocessed and trimmed, it is safe to set it to the raw length.

Shannon processes and gives assembled transcriptome at **output_path**. 
```
./Shannon_RNASeq_Cpp shannon -l 101 -s path_to_single_ended_read -o output_path -t 1 -g 0 -m 100G -u 4
```

For pair-ended fasta read file at **path_to_pair_ended_read_1** and **path_to_pair_ended_read_2**, with corresponding maximal **76** and **76** bases per read, Shannon processes and presents output at **output_path**. 
```
./Shannon_RNASeq_Cpp shannon -i 76 76 -p path_to_pair_ended_read_1 path_to_pair_ended_read_2 -o output_path -t 1 -g 0 -m 100G -u 4
```
**-o** When processing a fasta file for the first time, it is suggested to set the output_path to a non-existant new directory while working with command line. However, if a user wants to re-run the assembly with the same input file but different parameters, it is fine to use the same output directory, Shannon can automatically detect previously generated jellyfish file and sorted kmer file and load them. However, when kmer length is the changed parameter, the jellyfish file becomes invalid and hence a new output directory is needed.

**-t** specifies the number of running processes (1 threads per process). 

**-g** controls the read subsampling at components partition step. It has a great impact on the running speed when input read files are large, since the partitioned reads are later used for multibridging. Setting it to 0 disables subsample operation; otherwise a higher **g** tends to keep more reads. For read i, the sampling is based on the formula P_i(KEEP) = min(1, g/read_count_for_read_i), if it is not known, using 50-100 to avoid lossing too much reads (note it has same interpretation as trinity --normalize_max_read_cov inside in silico part)

**-m** specifies max amount of memory allowed for Shannon in the multi-bridge step, to prevent process crushing due to memory overflow. User can indicate the memory value using number of byte (like 1234), or using a number with unit (1k,1M,1G,1T). Notice, this memory limit is only effective in the multi-bridge step, and it manages the number of working processes to prevent memory crashing due to too many processes working at the same time.

**-u** number of parallel passed to linux sort function. Sorting a large file consumes a lot of memory, and linux kernel migth silently kill the process if the system is severely out of memory. Set the number to 8 if more than 100G is available. One way to check if the process is silently killed by system is to check if the file (kmer_unfiltered_sorted) is empty in the output directory.

### CMD\<partition>
Given a large noisy fasta-reads, CMD\<partition> performs error correction, and partitions the reads into multiple components. Each component contains two files, a kmer file which holds all error-corrected kmer and read files which selectively include original input reads for that component depending on the partition algorithm. 

```
./Shannon_RNASeq_Cpp partition -l read_length -s path_to_pair_ended_read_1 -o output_path
```
For all components, kmer files are stored together at **output_path/components_kmer**, read files are stored together at **output_path/components_reads**. 

### CMD\<multi-bridge>
Given one component (for example component 2) with kmer-reads pair as its input, it builds a de-bruijn graph using the kmer file. CMD\<multi-bridge> then condenses and pre-processes the graph by removing bubbles and suspicuous low abundance leaf nodes. After that it runs multi-bridge algorithm with the graph and the read files, which effectively reduces the number of xnodes within the graph if not all of them. 

```
./Shannon_RNASeq_Cpp multi-bridge -f path_to_kmer_file -l read_length -s path_to_read_file -o output_path 
```
Shannon uses three files for storing a connected graph:
```
node0   contains node in topological sorted order, where each node has a numerical ID
edge0   contains numerical ID pair to indicate a directed edge       
path0   a list of numerical ID to indicate that a read overlap such path       
where 0 is a tree component index
```
Here an index is introduced, since multi-bridge might break the graph into a forest. So for 10 tree components, 30 files are needed for storing the graph. Therefore For each graph component, three directories are used for storing such graph
```
output_path/comp_graph/Graph_edges_2/node(0,1,...,n)          
output_path/comp_graph/Graph_nodes_2/edge(0,1,...,n)                    
output_path/comp_graph/Graph_paths_2/path(0,1,...,n)          
where 2 is a graph component index
```

### CMD\<sparse-flow>
CMD\<sparse-flow> takes node,edge,path files as its inputs (for example node0,edge0,path0 ), then build a de-bruijn graph to run sparse flow algorithm. It generates reconstructed transcriptome.

```
./Shannon_RNASeq_Cpp sparse-flow --node_path path_to_node --edge_path path_to_edge --path_path path_to_path --output_path output_fasta_path
```

### CMD\<multi-graph>
CMD\<multi-graph> processes multiple components at the same time, use **-t** to specify the number of process where each process contains two threads, hence **-t 2** would allow 4 components running at the same time. Three modes are allowd
```
multi-bridge    run multi-bridge algorithm on multiple components
sparse-flow     run sparse-flow algorithm on all output graphs from multi-bridge
both            run multi-bridge and then sparse-flow on multiple components
```
The inputs to **multi-bridge** MODE are kmer-read inputs pair where kmer files for all components are stored together under **kmers_dir**, and all reads stored under **reads_dir**.

```
./Shannon_RNASeq_Cpp multi-graph multi-bridge -l single_read_length -c kmers_dir -r reads_dir -o output_path -t 1 -m 100G
```
The **sparse-flow** MODE, where **graph_dirs** holds all graphs
```
./Shannon_RNASeq_Cpp multi-graph sparse-flow -h graph_dirs -o output_path -t 1
```
The input to **both** mode is the same as **multi-bridge** MODE
```
./Shannon_RNASeq_Cpp multi-graph both -l single_read_length -c kmers_dir -r reads_dir -o output_path -t 1 -m 100G
```

### CMD\<find-rep>
CMD\<find-rep> reduces number of redundant outputs by checking for any two reconstructed contigs, if they has the similar length and both have the same starting r-mer and ending r-mer, where r is set to be 24.
```
./Shannon_RNASeq_Cpp find-rep --input path_to_input_fasta_file --output output_fasta_path
```
### CMD\<ref-align>
CMD\<ref-align> uses BLAT(https://genome.ucsc.edu/FAQ/FAQblat) to do the read alignment. It takes reconstructed transcriptome as its target, and for each reference transcript, it searches for the best match. To qualify for an alignment, the match length must be 90% of the reference length.

### CMD\<custom>
A developer interface, which takes a config file.

```
./Shannon_RNASeq_Cpp ref-align --input path_to_reconstructed_transcriptome --output_dir output_dir --ref path_to_reference_transcriptome
```

## Flowchart
The following flow chart illustrates how this implementation works:
* each cylinder represents a data storage, which is listed in the **output file structure** section. Item with [] means this file was just created by last functional block.
![Alt text](./doc/shannon_flowchart.png?raw=true "Optional Title")


## Output file structure
* algo_input (when only fasta file is provided, this directory takes jellyfish output)
	* kmer.dict
	* kmer.jf
* kmer_unfiltered_sorted
* kmer
* contig
* contig_filtered.fasta
* jf_stat
* components_kmer
	* kmer0, kmer1, kmer2 ...
* compomemts_read (# represents component index, & represents order of concatenation)
    * for single-ended reads
    	* comp0_0, comp0_1, ...
    	* comp1_0, comp1_1, ...
    	* ......
    	* compN_0, compN_1, ...
    * for paired-ended reads
    	* comp0_p1_0, comp0_p2_0, comp0_p1_1, comp0_p2_1 ... 
    	* ....
    	* compN_p1_0, compN_p2_0, compN_p1_1, compN_p2_1 ... 
* comp_array (a file for recording property of each component)
* comp_graph (# represents component index)
	* Graph_nodes_#
		* node0, node1, ... ,  single_node_#
	* Graph_edges_#
		* edge0, edge1, ...
	* Graph_paths_#
		* path0, path1, ...
* comp_graph_output (a temporary output file for storing reconstructed seq by each process)
* log_shannonC 
* reconstructed_seq.fasta_single (reconstructed seq from multibridge step)
* reconstructed_seq.fasta_sf (reconstructed seq from sparse flow step)
* reconstructed_seq.fasta (final output, combining filtered fasta_single, fasta_sf and contig_filtered.fasta)
* reconstructed_seq.fasta_unfiltered (simply combine fasta_single, fasta_sf and contig_filtered.fasta)
* mem_summary (a file summarizing memory usage)
* Shannon.timing ( a file keeping the time stamps once major tasks finish)


## Manual

### shannon
```
Usage: options_description [options]
Allowed options:
  --help                                produce help message
  -k [ --kmer_length ] arg (=25)        kmer length. default 25
  -d [ --double_strand ] arg (=1)       is reads considered double stranded.
  -l [ --single_read_length ] arg       single read length
  -i [ --pair_read_length ] arg         pair read length
  -s [ --SE_read_path ] arg             single ended read path 
  -p [ --PE_read_path ] arg             pair read path
  -e [ --output_seq_min_len ] arg (=200)
                                        minimal reconstructed output length
  -o [ --output_dir ] arg               output dir path
  --jf_path arg                         jf path from jellyfish
  --kmer_path arg                       kmer path from jellyfish
  --load_factor arg (=0.80000000000000004)
                                        kmer dictionary load factor, control 
                                        memory usage
  -u [ --num_sort_thread ] arg (=1)     argument passing to linux sort function
                                        option -parallel
  --sort_tmp_dir arg                    argument passing to linux sort function
                                        option -T which specifies a temporary
                                        directory used by linux sort function,
                                        /tmp wil be used if left empty
  --rmer_length arg (=15)               rmer length for error correction
  --threshold arg (=0.5)                ratio of common rmer between two 
                                        contigs for determing if two contifs 
                                        are repetative
  --min_count arg (=3)                  minimal kmer counts requirement for 
                                        serving as seeds in the greedy contig 
                                        formation step, also used in further 
                                        selection formula discussed in the 
                                        paper
  --min_length arg (=75)                minimal contigs length, also also used 
                                        in further selection formula discussed 
                                        in the paper
  --is_use_set arg (=0)                 use accummulating set for error 
                                        correction, it would reduce memory 
                                        usage, but give less performance
  --read_num_test arg (=3)              number test for a read to decide which 
                                        components the read should go to
  --is_assign_best arg (=1)             only assign a read to its best matching
                                        component or all components which work
  -g [ --read_sampler_k ] arg (=0)      given that many reads are redundant, 
                                        setting k to subsample matching reads 
                                        to speedup later works, where k can be 
                                        interpreted as coverage depth. k=0: no 
                                        subsample; k>0: probabilistic sample; 
                                        k<0: sample every (-k) read
  --metis_use_multiple_partition arg (=1)
                                        after one metis partition, reweight 
                                        graph so that previously broken edge 
                                        stays, it makes sure no information 
                                        loss due to partition, but generate 
                                        more output
  --metis_partition_size arg (=500)     metis parameter to control number of 
                                        contig within partitioned components
  --metis_non_partition_size arg (=500) for those graphs small enough not to be
                                        partitioned Shannon groups little 
                                        pieces together, and the parameter 
                                        control the collector component size
  --penalty arg (=5)                    metis parameter, for controlling 
                                        partition methods
  --overload arg (=2)                   metis parameter, for controlling 
                                        partition methods
  -r [ --reference ] arg                Reference path for the output to 
                                        compare with
  --rmer_length arg (=24)               rmer size used to check duplicates
  -m [ --avail_mem ] arg ()
                                        specify max memory for running
                                        multi-graph multi-bridge
  -t [ --num_process ] arg (=1)         Specify number processs for multi graph
                                        procedureswhere each process owns 
                                        one thread by default
  --max_hop_for_known_path arg (=30)    Max number of hops allowable for a read
                                        to cover
  --multiple_test arg (=8)              number of times that a linear program 
                                        sparse flow is solved

```

### partition
```
  --help                                produce help message
  -k [ --kmer_length ] arg (=25)        kmer length. default 25
  -d [ --double_strand ] arg (=1)       is reads considered double stranded.
  -l [ --single_read_length ] arg       single read length
  -i [ --pair_read_length ] arg         pair read length
  -s [ --SE_read_path ] arg             single ended read path 
  -p [ --PE_read_path ] arg             pair read path
  -e [ --output_seq_min_len ] arg (=200)
                                        minimal reconstructed output length
  -o [ --output_dir ] arg               output dir path
  --jf_path arg                         jf path from jellyfish
  --kmer_path arg                       kmer path from jellyfish
  --load_factor arg (=0.80000000000000004)
                                        kmer dictionary load factor, control 
                                        memory usage
  -u [ --num_sort_thread ] arg (=1)     argument passing to linux sort function
                                        option -parallel
  --sort_tmp_dir arg                    argument passing to linux sort function
                                        option -T which specifies a temporary
                                        directory used by linux sort function,
                                        /tmp wil be used if left empty
  --rmer_length arg (=15)               rmer length for error correction
  --threshold arg (=0.5)                ratio of common rmer between two 
                                        contigs for determing if two contifs 
                                        are repetative
  --min_count arg (=3)                  minimal kmer counts requirement for 
                                        serving as seeds in the greedy contig 
                                        formation step, also used in further 
                                        selection formula discussed in the 
                                        paper
  --min_length arg (=75)                minimal contigs length, also also used 
                                        in further selection formula discussed 
                                        in the paper
  --is_use_set arg (=0)                 use accummulating set for error 
                                        correction, it would reduce memory 
                                        usage, but give less performance
  --read_num_test arg (=3)              number test for a read to decide which 
                                        components the read should go to
  --is_assign_best arg (=1)             only assign a read to its best matching
                                        component or all components which work
  -g [ --read_sampler_k ] arg (=0)      given that many reads are redundant, 
                                        setting k to subsample matching reads 
                                        to speedup later works, where k can be 
                                        interpreted as coverage depth. k=0: no 
                                        subsample; k>0: probabilistic sample; 
                                        k<0: sample every (-k) read
  --metis_use_multiple_partition arg (=1)
                                        after one metis partition, reweight 
                                        graph so that previously broken edge 
                                        stays, it makes sure no information 
                                        loss due to partition, but generate 
                                        more output
  --metis_partition_size arg (=500)     metis parameter to control number of 
                                        contig within partitioned components
  --metis_non_partition_size arg (=500) for those graphs small enough not to be
                                        partitioned Shannon groups little 
                                        pieces together, and the parameter 
                                        control the collector component size
  --penalty arg (=5)                    metis parameter, for controlling 
                                        partition methods
  --overload arg (=2)                   metis parameter, for controlling 
                                        partition methods
```

### multi-bridge
```
Usage: options_description [options]
Allowed options:
  --help                                produce help message
  -f [ --kmer_path ] arg                kmer path used as trusted debruijn 
                                        graph edge
  -l [ --single_read_length ] arg       single read length
  -i [ --pair_read_length ] arg         pair read length
  -s [ --SE_read_path ] arg             single ended read path 
  -p [ --PE_read_path ] arg             pair read path
  --max_hop_for_known_path arg (=30)    Max number of hops allowable for a read
                                        to cover
  -e [ --output_seq_min_len ] arg (=200)
                                        minimal reconstructed output length
  -o [ --output_dir ] arg               output dir path
```

### sparse-flow
```
Usage: options_description [options]
Allowed options:
  --help                                produce help message
  --node_path arg                       input nodes path to sparse flow
  --edge_path arg                       input edges path to sparse flow
  --path_path arg                       input path path to sparse flow
  -e [ --output_seq_min_len ] arg (=200)
                                        minimal reconstructed output length
  --output_path arg                     output path
  --multiple_test arg (=8)              number of times that a linear program 
                                        sparse flow is solved

```

### multi-graph multi-bridge
```
Usage: options_description [options]
Allowed options:
  --help                                produce help message
  -k [ --kmer_length ] arg (=25)        kmer length. default 25
  -d [ --double_strand ] arg (=1)       is reads considered double stranded.
  -r [ --read_components_dir ] arg      specify reads_components dir produced 
                                        from partition step
  -c [ --kmer_components_dir ] arg      specify kmer_components dir produced 
                                        from partition step
  -m [ --avail_mem ] arg ()
                                        specify max memory for running
                                        multi-graph multi-bridge                            
  -l [ --single_read_length ] arg       single read length
  -i [ --pair_read_length ] arg         pair read length
  -e [ --output_seq_min_len ] arg (=200)
                                        minimal reconstructed output length
  -o [ --output_dir ] arg               output dir path
  --multiple_test arg (=8)              number of times that a linear program 
                                        sparse flow is solved
  --max_hop_for_known_path arg (=30)    Max number of hops allowable for a read
                                        to cover
  -t [ --num_process ] arg (=1)         Specify number processs for multi graph
                                        procedureswhere each process owns 
                                        one thread by default

```

### multi-graph sparse-flow
```
Usage: options_description [options]
Allowed options:
  --help                                produce help message
  -h [ --graphs_dir ] arg               specify graphs dir generated from 
                                        multi-bridge
  -e [ --output_seq_min_len ] arg (=200)
                                        minimal reconstructed output length
  -o [ --output_dir ] arg               output dir path
  --multiple_test arg (=8)              number of times that a linear program 
                                        sparse flow is solved
  -t [ --num_process ] arg (=1)         Specify number processs for multi graph
                                        procedureswhere each process owns 
                                        one thread by default
```

### find-rep 
```
Usage: options_description [options]
Allowed options:
  --help                                produce help message
  -e [ --output_seq_min_len ] arg (=200)
                                        minimal reconstructed output length
  --rmer_length arg (=24)               rmer size used to check duplicates
  --input arg                           input path
  --output arg                          output path
  -d [ --double_strand ] arg (=1)       Specify if reads are considered double 
                                        stranded.
```

### ref-align
```
Usage: options_description [options]
Allowed options:
  --help                produce help message
  --input arg           reconstucted fasta path
  --ref arg             reference path
  --output_dir arg      output dir

```



### custom
```
	* argument: a config file (More details in next section) 
	* custom setting provides a greater control over the whole processing, an user can start at any functional block to practice new parameters, provided that the corresponding inputs are available and placed at proper place [See Output file structure](#output-file-structure)
	* start from beginning, run the whole process: choose 12 (same as command line shannon)
	* run task A, B, C          	 : choose 3
	* run task B, C             	 : choose 14
	* run task B, C, D, E, F       	 : choose 21
	* run task C	            	 : choose 4	
	* run task C, D, E, F       	 : choose 20
	* run task D, E, F    	    	 : choose 13
	* run task D          	   	 : choose 9
	* run task E	    	         : choose 11
	* run single compoment in task D : choose 8, then enter component number
	* run single graph in task E	 : choose 10, then enter component component numner, then graph number
	* run task F                	 : choose 19	
	* run evaluation            	 : choose 16 (provided reference is given)
```

## Config file
Every task have various parameter, it is cumbumsome to have them all in argument. The config file serves to fit default value

Usage: each setting file should corresponds to a distinct output path

* kmer_setup
	* load_factor: load factor for the kmer dictionary created in task B
	* num_sort_thread: number sorting thread provided to linux sort function
	* rmer_length: an internal parameter for filtering noisy kmer (see paper)
	* threashold: an internal parameter for filtering noisy kmer (see paper)
	* min_count: an internal parameter for filtering noisy kmer (see paper)
	* min_length: an internal parameter for filtering noisy kmer (see paper)
	* is_use_set: usually set to false, instead use a dictionary filtering 

* contig_graph_setup
	* num_test: for each fasta-read (pair), it specifies number of kmer it would check for assigning the read to the components
	* is_assign_best: for num_test higher than 1, it determines if the read is assigned to one or all components, usually best is sufficient.
* metis_setup
	* partition_size: see metis manual
	* penalty: see metis manual
	* overload: see metis manual
	* use_multiple_partition: see paper
* multi_graph_setup
    * num_parallel: specify number of **process** created for multibridging all graphs. (each process would spawn two threads)
* seq_graph_setup
	* max_hop_for_known_path: known path is used in sparse flow, only consider path less than this number
* sparse_flow_setup
	* multiple_test: it specifies how many random trails are performed for solving L1 norm (see paper)
* find_rep_setup
    * rmer_length   

## Memory parameter
To address the memory issue, two several are defined 
* choose methods to represents kmer count: in the file src/shc_types.h
	* USE_APPROX_COUNT -- use 8 bytes for kmer dictionary value part
	* USE_EXACT_COUNT  -- use 16 bytes for kmer dictionary value part
	* experiment shows that for input with size 34G, USE_APPROX_COUNT takes 61G, USE_EXACT_COUNT takes 82G
* choose dictionary to store kmer count: in the file src/shc_setting.h
	* USE_HOPSCOTCH
	* USE_SPARSEPP
	* USE_HOPSCOTCH is faster in small size input, but takes too much memory if size is large and hence becomes slow due to internal memory management issue. USE_SPARSEPP is slightly slower, but is consistant




## Versioning
0.0.0

## Author

* **Bowen Xue** 

## License
Shannon is distributed under a GPL 3.0 license.

## Acknowledgments

