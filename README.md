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
* boost library (May be the most painful part)
	* 
* glpk 
	*  sudo apt-get install glpk

* hopscotch -- https://github.com/Tessil/hopscotch-map
* google sparsehash -- https://github.com/sparsehash/sparsehash
* sparsepp -- https://github.com/greg7mdp/sparsepp





## Running the tests

Two options: 1. shannon   2. custom
1. shannon inidicates a full process starting from breaking reads into kmers using jellyfish, to sparse flow reconstructing the final sequence
2. 19 detailed options which allow user to start from any steps

Usage: 
./Shannon_C_seq shannon options
./Shannon_C_seq custom setting_file


## Getting Started

* Try to run command **./Shannon_C_seq custom setting**



## Contributing


## Versioning


## Authors

* **Bowen Xue** 

## License


## Acknowledgments
