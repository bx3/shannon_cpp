# Shannon C++ RNA sequencing


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

C++ 11
boost 1.65.1 -- including non header portion -- http://www.boost.org/doc/libs/1_66_0/more/getting_started/unix-variants.html
metis 5.1.0 -- http://glaros.dtc.umn.edu/gkhome/metis/metis/download
hopscotch -- https://github.com/Tessil/hopscotch-map
google sparsehash -- https://github.com/sparsehash/sparsehash
sparsepp -- https://github.com/greg7mdp/sparsepp

glpk -- linux package
sort -- linux package

Tested with ubuntu 15.10

### Installing



## Running the tests

Two options: 1. shannon   2. custom
1. shannon inidicates a full process starting from breaking reads into kmers using jellyfish, to sparse flow reconstructing the final sequence
2. 19 detailed options which allow user to start from any steps

Usage: 
./Shannon_C_seq shannon options
./Shannon_C_seq custom setting_file


### And coding style tests



## Deployment



## Built With


## Contributing


## Versioning


## Authors

* **Bowen Xue** 

## License


## Acknowledgments
