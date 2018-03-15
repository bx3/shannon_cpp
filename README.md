# Shannon C++ RNA sequencing

This is a **C++** implementatin for the work **ShannonRNA: an Information-Optimal de Novo RNA-Seq Assembler**.  

A **python** implementation can be found at *http://sreeramkannan.github.io/Shannon/*, which served as a foundation for the current **C++** implementation.  

However, this C++ version significantly improves both time and memory efficiency, roughly 10x, so you would prefer using this one.


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
* boost library (May be the most painful part if you have an old one on system)
	* follow instruction from [boost](http://www.boost.org/doc/libs/1_66_0/more/getting_started/unix-variants.html)
	* remember to install the non-header part of boost. This implementation needs
		* Boost.ProgramOptions
		* Boost.Filesystem 
* glpk 
	*  sudo apt-get install libglpk-dev
	*  OR follow instructions from [GLPK](https://en.wikibooks.org/wiki/GLPK/Linux_OS)
* google sparsehash 
	* follow instructions from [google sparsehash](https://github.com/sparsehash/sparsehash)
* hopscotch -- https://github.com/Tessil/hopscotch-map
* sparsepp -- https://github.com/greg7mdp/sparsepp

## Getting Started

* Try to run command **./Shannon_C_seq custom setting**


## Output Files structure

./  
1. algo_input (when only fasta file is provided, this folder becomes jellyfish output)
	* kmer.dict
	* kmer.jf
2. kmer_unfiltered_sorted
3. kmer
4. contig
5. jf_stat
6. components_kmer
	* kmer# (is partitioned kmer file for #-th component)	
7. compomemts_read
	* comp#_@ (is partititoned fasts files for #-th components, @ is the chained filed number)
	* comp#_p1_@ (used for paired read )
	* comp#_p2_@ (used for paired read )
8. comp_array (a file for recording property of each component)
9. comp_graph
	* Graph_nodes_#
		* node@ 
		* single_node_#
	* Graph_edges_#
		* edge@
	* Graph_paths_E
		* path@
10. comp_graph_output (a temporary output file for storing reconstructed seq by each process)
11. log_shannonC 
12. reconstructed_seq.fasta (final output, combine fasta_single and fasta_sf, then perform filtering)
13. reconstructed_seq.fasta_single (reconstructed seq from multibridge step)
14. reconstructed_seq.fasta_sf (reconstructed seq from sparse flow step)
15. reconstructed_seq.fasta_unfiltered (simply combine fasta_single and fasta_sf)

## How things work
The shannon consists of four major components:
1. noisy kmer filtering
2. contig graph partitioning
3. sequence graph multibridging
4. sparse flowing path decomposing

![Alt text](./doc/?raw=true "Optional Title")



## Command line explain
Two modes are available 
* shannon
	* shannon provide a way to specify most important parameters through command line 
 	* required field
 		* j -- config file, specifying all detail parameter (More details in next section)
 		* k -- kmer length
 		* o -- output directory
 		* input fasta path and lenght
 			* s -- single ended fasta path
 			* l -- single ended fasta read typical length
 			* p -- paired ended fasta path
 			* i -- paired ended fasta read typical length in a list
 			* or include both
 	* optional field
 		* it is fine only provide fasta path, but if jf and kmer dictionary is available, can directly use them
 			* n -- jf path
 			* m -- kmer path
 		* if reference is available for evaluation
 			* r -- reference 
 		* d -- is single stranded
* custom 
	* custom provide a fine-tuned control over the whole Shannon_RNA_Seq processing, user can start at any check point to practice new parameters setting, provided that input are available
	* 19 tests are available
	* 





## Running the tests






## Contributing


## Versioning


## Authors

* **Bowen Xue** 

## License


## Acknowledgments
