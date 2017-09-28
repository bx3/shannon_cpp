/* 
 * File:   shc_google_sparsehash.h
 * Author: bx
 *
 * Created on September 3, 2017, 2:48 PM
 */

#ifndef SHC_GOOGLE_SPARSEHASH_H
#define	SHC_GOOGLE_SPARSEHASH_H
#include <stdint.h>
#include "shc_type.h"
#include <google/sparse_hash_map>
#include <google/dense_hash_map>
#include <google/dense_hash_set>
//using sparse_hash_map; 

struct equ64 {

  bool operator()(const uint64_t& kmer_val_a, const uint64_t& kmer_val_b) const {
	return( kmer_val_a == kmer_val_b);
  }
};


struct hash_u64 {
  
  uint64_t operator() (const uint64_t& kmer_val) const {

	return(kmer_val);
  }
};

struct eq_contig {

  bool operator()(const contig_num_t& kmer_val_a, const contig_num_t& kmer_val_b) const {
	return( kmer_val_a == kmer_val_b);
  }
};

struct hash_contig {
  
  uint64_t operator() (const contig_num_t& kmer_val) const {

	return(kmer_val);
  }
};


typedef google::dense_hash_map<uint64_t, Kmer_info, hash_u64, equ64> Kmer_counter_map;
typedef google::dense_hash_map<uint64_t, Kmer_info, hash_u64, equ64>::iterator Kmer_counter_map_iterator;
typedef google::dense_hash_map<uint64_t, Kmer_info, hash_u64, equ64>::const_iterator Kmer_counter_map_const_iterator;

typedef google::sparse_hash_map<uint64_t, std::vector<contig_num_t>, hash_u64, equ64> K1mer_contig_map;
typedef google::sparse_hash_map<uint64_t, std::vector<contig_num_t>, hash_u64, equ64>::iterator K1mer_contig_map_iterator;
typedef google::sparse_hash_map<uint64_t, std::vector<contig_num_t>, hash_u64, equ64>::const_iterator K1mer_contig_map_const_iterator;

typedef google::sparse_hash_map<uint64_t, std::vector<contig_num_t>, hash_u64, equ64> Rmer_contig_map;
typedef google::sparse_hash_map<uint64_t, std::vector<contig_num_t>, hash_u64, equ64>::iterator Rmer_contig_map_iterator;
typedef google::sparse_hash_map<uint64_t, std::vector<contig_num_t>, hash_u64, equ64>::const_iterator Rmer_contig_map_const_iterator;

typedef google::sparse_hash_map<uint64_t, contig_num_t, hash_u64, equ64> Rmer_count_map;
typedef google::sparse_hash_map<uint64_t, contig_num_t, hash_u64, equ64>::iterator Rmer_count_map_iterator;
typedef google::sparse_hash_map<uint64_t, contig_num_t, hash_u64, equ64>::const_iterator Rmer_count_map_const_iterator;

typedef google::sparse_hash_map<contig_num_t, rmer_count_t, hash_contig, eq_contig> Contig_count_map;
typedef google::sparse_hash_map<contig_num_t, rmer_count_t, hash_contig, eq_contig>::iterator Contig_count_map_iterator;
typedef google::sparse_hash_map<contig_num_t, rmer_count_t, hash_contig, eq_contig>::const_iterator Contig_count_map_const_iterator;

//typedef google::sparse_hash_map<uint64_t, contig_num_t, hashme_K1, eqstr_K1> Rmer_counter_map;
typedef google::dense_hash_set<contig_num_t, hash_contig, eq_contig> Contig_set;
typedef google::dense_hash_set<contig_num_t, hash_contig, eq_contig>::iterator Contig_set_iterator;
typedef google::dense_hash_set<contig_num_t, hash_contig, eq_contig>::const_iterator Contig_set_const_iterator;

#endif	/* SHC_GOOGLE_SPARSEHASH_H */

