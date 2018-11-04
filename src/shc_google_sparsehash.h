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

#include "hopscotch_map.h"
#include "hopscotch_set.h"

#include <sparsepp/spp.h>
#include "encoding.h"
#include "shc_setting.h"


struct equ64 {

  bool operator()(const uint64_t& kmer_val_a, const uint64_t& kmer_val_b) const {
	return( kmer_val_a == kmer_val_b);
  }
};

struct hash_dna_seq {
    uint64_t operator() (const std::string& seq) const {
       //return (u ^ 14695981039346656037ULL) * 1099511628211ULL;
       uint64_t n = 1;
       for(int i=0; i<std::ceil(seq.size()/32); i++)
       {
            uint64_t byte;
            if(i!=std::ceil(seq.size()/32)-1)
                encode_kmer(&seq.at(i*32), &byte, 32);
            else
                encode_kmer(&seq.at(i*32), &byte, seq.size()%32);
            n *= byte;
       }


      return n;

      //size_t operator()(uint64_t k) const { return (k ^ 14695981039346656037ULL) * 1099511628211ULL; }
    }
};


struct hash_u64 {

  uint64_t operator() (const uint64_t& u) const {
     //return (u ^ 14695981039346656037ULL) * 1099511628211ULL;
    //const uin32_t base = 2166136261U;
    //const uin32_t p = 16777619UL;

    //size_t res = 2166136261U;
    //res ^= u;
    //res *= 16777619UL;


    uint64_t v = u * 3935559000370003845 + 2691343689449507681;

    v ^= v >> 21;
    v ^= v << 37;
    v ^= v >>  4;

    v *= 4768777513237032717;

    v ^= v << 20;
    v ^= v >> 41;
    v ^= v <<  5;

    return v;

    //size_t operator()(uint64_t k) const { return (k ^ 14695981039346656037ULL) * 1099511628211ULL; }
  }
};

struct equ32 {

  bool operator()(const uint32_t& kmer_val_a, const uint32_t& kmer_val_b) const {
	return( kmer_val_a == kmer_val_b);
  }
};


struct hash_u32 {

  uint32_t operator() (const uint32_t& hash) const {
      //size_t operator()(uint32_t k) const { return (k ^ 2166136261U)  * 16777619UL; }
      uint32_t v = (hash ^ 2166136261U)  * 16777619UL;
	v += v << 3;
        v ^= v >> 11;
        v += v << 15;

      return(v);
  }
};

struct eq_contig {

  bool operator()(const contig_num_t& kmer_val_a, const contig_num_t& kmer_val_b) const {
	return( kmer_val_a == kmer_val_b);
  }
};

struct hash_contig {

    uint64_t operator() (const uint64_t& u) const {
       //return (u ^ 14695981039346656037ULL) * 1099511628211ULL;

      uint64_t v = u * 3935559000370003845 + 2691343689449507681;

      v ^= v >> 21;
      v ^= v << 37;
      v ^= v >>  4;

      v *= 4768777513237032717;

      v ^= v << 20;
      v ^= v >> 41;
      v ^= v <<  5;

      return v;
    //size_t operator()(uint64_t k) const { return (k ^ 14695981039346656037ULL) * 1099511628211ULL; }
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

#ifdef USE_DENSE_KMER
typedef google::dense_hash_map<uint64_t, Kmer_info, hash_u64, equ64> Kmer_counter_map;
typedef google::dense_hash_map<uint64_t, Kmer_info, hash_u64, equ64>::iterator Kmer_counter_map_iterator;
typedef google::dense_hash_map<uint64_t, Kmer_info, hash_u64, equ64>::const_iterator Kmer_counter_map_const_iterator;

#define RESERVE_NUM_KMER(dict, size)           dict.resize(size);
#define ITER_GET_VALUE(iter)                   iter->second
#define ITER_SET_CONTIG(iter, contig_)         iter->second.contig = contig_;
#define ITER_SET_INFO(iter, info_)             iter->second.info = info_;
#define ITER_SET_USED(iter, used_)             iter->second.used = used_;
#define INIT_DICT(dict, empty_key, delete_key) \
        do {dict.set_deleted_key(delete_key); dict.set_empty_key(empty_key);} while (0)//also work for set

typedef google::dense_hash_set<uint64_t, hash_u64, equ64> Rmer_set;
typedef google::dense_hash_set<uint64_t, hash_u64, equ64>::iterator Rmer_set_iterator;
typedef google::dense_hash_set<uint64_t, hash_u64, equ64>::const_iterator Rmer_set_const_iterator;

typedef google::dense_hash_map<uint64_t, std::vector<contig_num_t>, hash_u64, equ64> Rmer_contig_map;
typedef google::dense_hash_map<uint64_t, std::vector<contig_num_t>, hash_u64, equ64>::iterator Rmer_contig_map_iterator;
typedef google::dense_hash_map<uint64_t, std::vector<contig_num_t>, hash_u64, equ64>::const_iterator Rmer_contig_map_const_iterator;

#endif

#ifdef USE_SPARSE_KMER
typedef google::sparse_hash_map<uint64_t, Kmer_info, hash_u64, equ64> Kmer_counter_map;
typedef google::sparse_hash_map<uint64_t, Kmer_info, hash_u64, equ64>::iterator Kmer_counter_map_iterator;
typedef google::sparse_hash_map<uint64_t, Kmer_info, hash_u64, equ64>::const_iterator Kmer_counter_map_const_iterator;


#define RESERVE_NUM_KMER(dict, size)          dict.resize(size);
#define ITER_GET_VALUE(iter)                  iter->second
#define ITER_SET_CONTIG(iter, contig_)        iter->second.contig = contig_;
#define ITER_SET_INFO(iter, info_)            iter->second.info = info_;
#define ITER_SET_USED(iter, used_)            iter->second.used = used_;
#define INIT_DICT(dict, empty_key, delete_key) \
        do { dict.set_deleted_key(delete_key);} while(0)

typedef google::dense_hash_set<uint64_t, hash_u64, equ64> Rmer_set;
typedef google::dense_hash_set<uint64_t, hash_u64, equ64>::iterator Rmer_set_iterator;
typedef google::dense_hash_set<uint64_t, hash_u64, equ64>::const_iterator Rmer_set_const_iterator;
#endif

#ifdef USE_HOPSCOTCH
typedef tsl::hopscotch_map<uint64_t, Kmer_info, hash_u64, equ64> Kmer_counter_map;
typedef tsl::hopscotch_map<uint64_t, Kmer_info, hash_u64, equ64>::iterator Kmer_counter_map_iterator;
typedef tsl::hopscotch_map<uint64_t, Kmer_info, hash_u64, equ64>::const_iterator Kmer_counter_map_const_iterator;

#define RESERVE_NUM_KMER(dict, size)                dict.reserve(size);
#define ITER_INCREMENT_COUNT(iter, count_)    iter.value().count = ADD_COUNT((iter.value().count) , count_);
#define ITER_GET_VALUE(iter)                  iter.value()
#define ITER_SET_CONTIG(iter, contig_)        iter.value().contig = contig_;
#define ITER_SET_INFO(iter, info_)            iter.value().info = info_;
#define ITER_SET_USED(iter, used_)                   iter.value().used = used_;
#define INIT_DICT(dict, empty_key, delete_key) \
        do {;} while(0)

typedef tsl::hopscotch_set<uint64_t, hash_u64, equ64> Rmer_set;
typedef tsl::hopscotch_set<uint64_t, hash_u64, equ64>::iterator Rmer_set_iterator;
typedef tsl::hopscotch_set<uint64_t, hash_u64, equ64>::const_iterator Rmer_set_const_iterator;

typedef tsl::hopscotch_map<uint64_t, std::vector<contig_num_t>, hash_u64, equ64> Rmer_contig_map;
typedef tsl::hopscotch_map<uint64_t, std::vector<contig_num_t>, hash_u64, equ64>::iterator Rmer_contig_map_iterator;
typedef tsl::hopscotch_map<uint64_t, std::vector<contig_num_t>, hash_u64, equ64>::const_iterator Rmer_contig_map_const_iterator;
#endif

#ifdef USE_SPARSEPP
typedef spp::sparse_hash_map<uint64_t, Kmer_info, hash_u64, equ64> Kmer_counter_map;
typedef spp::sparse_hash_map<uint64_t, Kmer_info, hash_u64, equ64>::iterator Kmer_counter_map_iterator;
typedef spp::sparse_hash_map<uint64_t, Kmer_info, hash_u64, equ64>::const_iterator Kmer_counter_map_const_iterator;

#define RESERVE_NUM_KMER(dict, size)                dict.reserve(size);
#define ITER_GET_VALUE(iter)                  iter->second
#define ITER_INCREMENT_COUNT(iter, count_)    iter->second.count = ADD_COUNT(iter->second.count ,count_);
#define ITER_SET_CONTIG(iter, contig_)        iter->second.contig = contig_;
#define ITER_SET_INFO(iter, info_)            iter->second.info = info_;
#define ITER_SET_USED(iter, used_)            iter->second.used = used_;
#define INIT_DICT(dict, empty_key, delete_key) \
        do {dict.set_deleted_key(delete_key);} while (0)

typedef tsl::hopscotch_set<uint64_t, hash_u64, equ64> Rmer_set;
typedef tsl::hopscotch_set<uint64_t, hash_u64, equ64>::iterator Rmer_set_iterator;
typedef tsl::hopscotch_set<uint64_t, hash_u64, equ64>::const_iterator Rmer_set_const_iterator;

typedef spp::sparse_hash_map<uint64_t, std::vector<contig_num_t>, hash_u64, equ64> Rmer_contig_map;
typedef spp::sparse_hash_map<uint64_t, std::vector<contig_num_t>, hash_u64, equ64>::iterator Rmer_contig_map_iterator;
typedef spp::sparse_hash_map<uint64_t, std::vector<contig_num_t>, hash_u64, equ64>::const_iterator Rmer_contig_map_const_iterator;

#endif


typedef google::sparse_hash_map<uint64_t, std::vector<contig_num_t>, hash_u64, equ64> K1mer_contig_map;
typedef google::sparse_hash_map<uint64_t, std::vector<contig_num_t>, hash_u64, equ64>::iterator K1mer_contig_map_iterator;
typedef google::sparse_hash_map<uint64_t, std::vector<contig_num_t>, hash_u64, equ64>::const_iterator K1mer_contig_map_const_iterator;

typedef google::sparse_hash_map<uint64_t, contig_num_t, hash_u64, equ64> Rmer_count_map;
typedef google::sparse_hash_map<uint64_t, contig_num_t, hash_u64, equ64>::iterator Rmer_count_map_iterator;
typedef google::sparse_hash_map<uint64_t, contig_num_t, hash_u64, equ64>::const_iterator Rmer_count_map_const_iterator;



typedef google::dense_hash_map<contig_num_t, rmer_count_t, hash_contig, eq_contig> Contig_count_map;
typedef google::dense_hash_map<contig_num_t, rmer_count_t, hash_contig, eq_contig>::iterator Contig_count_map_iterator;
typedef google::dense_hash_map<contig_num_t, rmer_count_t, hash_contig, eq_contig>::const_iterator Contig_count_map_const_iterator;

//typedef google::sparse_hash_map<uint64_t, contig_num_t, hashme_K1, eqstr_K1> Rmer_counter_map;
typedef google::dense_hash_set<contig_num_t, hash_contig, eq_contig> Contig_set;
typedef google::dense_hash_set<contig_num_t, hash_contig, eq_contig>::iterator Contig_set_iterator;
typedef google::dense_hash_set<contig_num_t, hash_contig, eq_contig>::const_iterator Contig_set_const_iterator;

typedef google::dense_hash_map<read_num_t, read_length_t, hash_read, eq_read> Read_len_map;
typedef google::dense_hash_map<read_num_t, read_length_t, hash_read, eq_read>::iterator Read_len_map_iterator;
typedef google::dense_hash_map<read_num_t, read_length_t, hash_read, eq_read>::const_iterator Read_len_map_const_iterator;
#endif	/* SHC_GOOGLE_SPARSEHASH_H */
