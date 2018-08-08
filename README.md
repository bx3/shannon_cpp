# Shannon C++ RNA sequencing

This is a **C++** implementatin for the work **ShannonRNA: an Information-Optimal de Novo RNA-Seq Assembler**.  A **python** implementation can be found at *http://sreeramkannan.github.io/Shannon/*, which served as a foundation for the current **C++** implementation. However, this C++ version significantly improves both time and memory efficiency, roughly 10x, so you would prefer using this one.

## Content
- [Prerequisites](#prerequisites)
- [Installing](#installing)
- [Input requirement](#input-requirement)
- [Getting started](#getting-started)
- [Output file structure](#output-file-structure)
- [Flowchart](#flowchart)
- [Usage](#usage)
- [Config file](#config-file)
- [Memory parameter](#memory-parameter)
- [Versioning](#versioning)
- [Author](#author)
- [License](#license)
- [Acknowledgments](#acknowledgments)
### Prerequisites

* ubuntu 15.10 or higher
* g++

### Installing
* Install METIS
	* sudo apt-get install libmetis-dev
	* OR follow instructions from [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/download) .
* Install cmake 
	* sudo apt-get install cmake
* Install Jellyfish
	* follow instrctions from [JellyFish](http://www.cbcb.umd.edu/software/jellyfish/) 
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
### CMD<shannon>
Shannon performs denovo transcriptome assembly by first partitioning reads into components, which can be individually run with CMD<partition>. Shannon then builds a de-bruijn graph for each component and runs CMD<multi-bridge> and CMD<sparse-flow> on that graph. To process all components at the same time, CMD<multi-graph> spawns multiple processes to facilitate parallel computing. At this step, Shannon has obtained reconstructed transcriptome and it uses CMD<find-rep> to remove redundancy. After all, CMD<ref-align> is optional, and allows comparison to the *ground true* result with BLAT. 
The Input to CMD<shannon> takes a large number of raw reads (or Rcorrected reads), and Shannon runs those subcommands and store both intermidiate results and final output at output directory.

For a single-ended fasta read file at **path_to_single_ended_read**, with accepted 101 bases per read, Shannon processes and gives output at **output_path** 
```
./Shannon_RNASeq_Cpp shannon -l 101 -s path_to_single_ended_read -o output_path
```

For pair-ended fasta read file at **path_to_pair_ended_read_1** and **path_to_pair_ended_read_2**, with corresponding accepted **76** and **76** bases per read, Shannon processes and gives output at **output_path**
```
./Shannon_RNASeq_Cpp shannon -i 76 76 -p path_to_pair_ended_read_1 path_to_pair_ended_read_2 -o output_path
```

### CMD<partition>
Given a large noisy fasta-reads, CMD<partition> performs error correction, and partitions the reads into multiple components. Each component contains two files, a kmer file which holds all error-corrected kmer and read files which selectively include original input reads for that component depending on the partition algorithm. 

```
./Shannon_RNASeq_Cpp partition -l read_length -s path_to_pair_ended_read_1 -o output_path
```
For all components, kmer files are stored together at **output_path/components_kmer**, read files are stored together at **output_path/components_reads**. 

### CMD<multi-bridge>
Given one component (for example component 2) with kmer-reads pair as its input, it builds a de-bruijn graph using the kmer file. CMD<multi-bridge> then condenses and pre-processes the graph by removing bubbles and suspicuous low abundance leaf nodes. After that it runs multi-bridge algorithm with the graph and the read files, which effectively reduces the number of xnodes within the graph if not all of them. 

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

### CMD<sparse-flow>
CMD<sparse-flow> takes node,edge,path files as its inputs (for example node0,edge0,path0 ), then build a de-bruijn graph to run sparse flow algorithm. It generates reconstructed transcriptome.

```
./Shannon_RNASeq_Cpp sparse-flow --node_path path_to_node --edge_path path_to_edge --path_path path_to_path --output_path output_fasta_path
```

### CMD<multi-graph>
CMD<multi-graph> processes multiple components at the same time, use **-t** to specify the number of process where each process contains two threads, hence **-t 2** would allow 4 components running at the same time. Three modes are allowd
```
multi-bridge    run multi-bridge algorithm on multiple components
sparse-flow     run sparse-flow algorithm on all output graphs from multi-bridge
both            run multi-bridge and then sparse-flow on multiple components
```
The inputs to **multi-bridge** MODE are kmer-read inputs pair where kmer files for all components are stored together under **kmers_dir**, and all reads stored under **reads_dir**.

```
./Shannon_RNASeq_Cpp multi-graph multi-bridge -l single_read_length -c kmers_dir -r reads_dir -o output_path -t 1
```
The **sparse-flow** MODE, where **graph_dirs** holds all graphs
```
./Shannon_RNASeq_Cpp multi-graph sparse-flow -h graph_dirs -o output_path -t 1
```
The input to **both** mode is the same as **multi-bridge** MODE
```
./Shannon_RNASeq_Cpp multi-graph both -l single_read_length -c kmers_dir -r reads_dir -o output_path -t 1
```

### CMD<find-rep>
CMD<find-rep> reduces number of redundant outputs by checking for any two reconstructed contigs, if they has the similar length and both have the same starting r-mer and ending r-mer, where r is set to be 24.
```
./Shannon_RNASeq_Cpp find-rep --input path_to_input_fasta_file --output output_fasta_path
```
### CMD<ref-align>
CMD<ref-align> uses BLAT(https://genome.ucsc.edu/FAQ/FAQblat) to do the read alignment. It takes reconstructed transcriptome as its target, and for each reference transcript, it searches for the best match. To qualify for an alignment, the match length must be 90% of the reference length.

```
./Shannon_RNASeq_Cpp ref-align --input path_to_reconstructed_transcriptome --output_dir output_dir --ref path_to_reference_transcriptome
```

## Output file structure

./  
* algo_input (when only fasta file is provided, this directory takes jellyfish output)
	* kmer.dict
	* kmer.jf
* kmer_unfiltered_sorted
* kmer
* contig
* jf_stat
* components_kmer
	* kmer# (is partitioned kmer file for #-th component)	
* compomemts_read (# represents component index, & represents order of concatenation)
	* comp#_& (is partititoned fasts files for #-th components, @ is the chained filed number)
	* comp#_p1_& (used for paired read )
	* comp#_p2_& (used for paired read )
* comp_array (a file for recording property of each component)
* comp_graph (# represents component index, @ represents graph index)
	* Graph_nodes_#
		* node@ 
		* single_node_#
	* Graph_edges_#
		* edge@
	* Graph_paths_#
		* path@
* comp_graph_output (a temporary output file for storing reconstructed seq by each process)
* log_shannonC 
* reconstructed_seq.fasta (final output, combine fasta_single and fasta_sf, then perform filtering)
* reconstructed_seq.fasta_single (reconstructed seq from multibridge step)
* reconstructed_seq.fasta_sf (reconstructed seq from sparse flow step)
* reconstructed_seq.fasta_unfiltered (simply combine fasta_single and fasta_sf)

## Flowchart
The following flow chart illustrates how this implementation works:
* each cylinder represents inputs, which was listed in the **output file structure** above, required by each functional block. Item with [] means this file was just created by last functional block.
![Alt text](./doc/shannon_flowchart.PNG?raw=true "Optional Title")


## Usage
Two modes are available 
* shannon
	* shannon provide a way to specify most important parameters through command line 
 	* required field
 		* j -- config file, specifying all detail parameter (More details in next section)
 		* k -- kmer length (**<=32**)
 		* o -- output directory
 		* input fasta path and lenght
			* single ended input
 				* s -- single ended fasta path
 				* l -- single ended fasta read typical length
			* paired ended input
 				* p -- paired ended fasta path
 				* i -- paired ended fasta read typical length in a list
 			* OR include both
 	* optional field
 		* it is fine only provide fasta path, but if jf and kmer dictionary is available, and user does not want to do duplicate work, paths to those file can be provided by
 			* n -- jf path
 			* m -- kmer path
 		* if reference is available for evaluation
 			* r -- reference 
 		* d -- is single stranded
* custom 
	* argument: a config file (More details in next section) 
	* custom setting provides a fine-tuned control over the whole processing, an user can start at any functional block to practice new parameters
	* 19 tests are available to start at any point in the flowchart, provided that the corresponding inputs inside a cylinder are available and placed at proper place [See Output file structure](#output-file-structure)
	* start from beginning, run the whole process: choose 12 (same as command line)
	* run task A, B, C          	 : choose 3
	* run task B, C             	 : choose 14
	* run task C	            	 : choose 4	
	* run task C, D, E, F       	 : choose 20
	* run task D, E, F    	    	 : choose 13
	* run task D          	   	 : choose 9
	* run task E	    	         : choose 11
	* run single compoment in task D : choose 8, then enter component number
	* run single graph in task E	 : choose 10, then enter component component numner, then graph number
	* run task F                	 : choose 19	
	* run evaluation            	 : choose 16 (provided reference is given)


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
	* is_use_set: usually set to false, instead use a dictionary filtering (see paper)

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
	* max_hop_for_pair_search: given a pair read, it would search for a unique path(if available) between the paired end and constructed a pseudo super-read
		* negative number: turn it off
		* 0: consider only if two termainl ends of the paired read both reside on the same node
		* positive number: check for how many hops
	* mate_pair_len: for the paired search mentioned above, only accept the super-read if middle part has a length less than this number
	* single_node_seq_len: when a graph is too well multibridged, it would yield many single nodes, consider this node to be a reconstructed RNA sequence iff its length is higher than this number
	* max_hop_for_known_path: known path is used in sparse flow, only consider path less than this number
* sparse_flow_setup
	* multiple_test: it specifies how many random trails are performed for solving L1 norm (see paper)
	* sf_num_parallel: number of **process** created for solving sparse flow (each process implicitly creats two threads)
	

## Memory parameter
To address the memory issue, two several are defined 
* choose methods to represents kmer count: in the file src/shc_types.h
	* USE_APPROX_COUNT -- use 8 bytes for kmer dictionary value part
	* USE_EXACT_COUNT  -- use 16 bytes for kmer dictionary value part
	* experiment shows thatn for input with size 34G, USE_APPROX_COUNT takes 61G, USE_EXACT_COUNT takes 82G
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

