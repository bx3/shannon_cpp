/*
 * File:   shc_google_sparsehash.h
 * Author: bx
 *
 * Created on September 3, 2017, 2:48 PM
 */

#ifndef SHC_GOOGLE_SPARSEHASH_H
#define	SHC_GOOGLE_SPARSEHASH_H
#include <stdint.h>
#include <unordered_map>
#include "shc_type.h"
#include <google/sparse_hash_map>
#include <google/sparse_hash_set>
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

struct equ32 {

  bool operator()(const uint32_t& kmer_val_a, const uint32_t& kmer_val_b) const {
	return( kmer_val_a == kmer_val_b);
  }
};


struct hash_u32 {

  uint64_t operator() (const uint32_t& kmer_val) const {

	return(kmer_val);
  }
};

struct eq_contig {

  bool operator()(const contig_num_t& kmer_val_a, const contig_num_t& kmer_val_b) const {
	return( kmer_val_a == kmer_val_b);
  }
};

struct hash_contig {

  contig_num_t operator() (const contig_num_t& kmer_val) const {

	return(kmer_val);
  }
};

struct eq_read {

  bool operator()(const read_num_t& kmer_val_a, const read_num_t& kmer_val_b) const {
	return( kmer_val_a == kmer_val_b);
  }
};

struct hash_read {

  read_num_t operator() (const read_num_t& kmer_val) const {

	return(kmer_val);
  }
};

struct eqstr
{
  bool operator()(const std::string s1, const std::string s2) const
  {
    return (s1 == s2);
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

typedef google::dense_hash_map<read_num_t, read_length_t, hash_read, eq_read> Read_len_map;
typedef google::dense_hash_map<read_num_t, read_length_t, hash_read, eq_read>::iterator Read_len_map_iterator;
typedef google::dense_hash_map<read_num_t, read_length_t, hash_read, eq_read>::const_iterator Read_len_map_const_iterator;
#endif	/* SHC_GOOGLE_SPARSEHASH_H */
