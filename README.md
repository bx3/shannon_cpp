# Shannon C++ RNA sequencing

This is a **C++** implementatin for the work **ShannonRNA: an Information-Optimal de Novo RNA-Seq Assembler**.  

A **python** implementation can be found at *http://sreeramkannan.github.io/Shannon/*, which served as a foundation for the current **C++** implementation.  

However, this C++ version significantly improves both time and memory efficiency, roughly 10x, so you would prefer using this one.

## Content
- [Prerequisites](#prerequisites)
- [Installing](#installing)
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
* C++ 11

### Installing

* Install METIS
	* sudo apt-get install libmetis-dev
	* OR follow instructions from [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/download) .
* Install cmake 
	* sudo apt-get install cmake
* Install Jellyfish
	* follow instrctions from [JellyFish](http://www.cbcb.umd.edu/software/jellyfish/) 
* boost library (May be the most painful part if you have an old one on system, cmake might not link properly)
	* follow instruction from [boost](http://www.boost.org/doc/libs/1_66_0/more/getting_started/unix-variants.html)
	* ~~sudo apt-get install libboost-all-dev~~, since it would install an older 
	
	, which would not compile
	* remember to install the non-header part of boost. This implementation needs
		* Boost.ProgramOptions
		* Boost.Filesystem 
		* Boost.System
* glpk 	
	*  install version 4.62 or higher, follow instructions from [GLPK](https://en.wikibooks.org/wiki/GLPK/Linux_OS)
	*  ~~sudo apt-get install libglpk-dev~~, since it would install an older version, which gives error
* google sparsehash 
	* follow instructions from [google sparsehash](https://github.com/sparsehash/sparsehash)
* hopscotch -- https://github.com/Tessil/hopscotch-map (already included in the repo, no action needed)
* sparsepp -- https://github.com/greg7mdp/sparsepp
	* download from link above, extract and place the whole files in the cmake directory

## Getting Started
* Example, there is a single ended read (at dir single_read_path) with read length 100, save it to ouput_path to 
* Try to run command with absolute paths
* **./Shannon_C_seq shannon -l default_setting -k 25 -l 100 -s single_read_path -o output_path**

## Output file structure

./  
* algo_input (when only fasta file is provided, this folder becomes jellyfish output)
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
 		* k -- kmer length
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
	* start from beginning, run the whole process: choose 12
	* run task A, B, C: choose 3
	* run task B, C   : choose 14
	* run task C	  : choose 4	
	* run task D, E, F: choose 13
	* run task D      : choose 9
	* run task E	  : choose 11
	* run single graph in task D: choose 8, then enter graph number
	* run single graph in task E: choose 10, then enter component graph numner, then graph number
	* run task F      : choose 19
	* run evaluation  : choose 16 (provided reference is given)


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
1.0.0

## Author

* **Bowen Xue** 

## License
Shannon is distributed under a GPL 3.0 license.

## Acknowledgments
