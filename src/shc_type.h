/*
 * File:   shc_type.h
 * Author: bx
 *
 * Created on September 2, 2017, 11:54 PM
 */

#ifndef SHC_TYPE_H
#define	SHC_TYPE_H

#include <stdint.h>
#include <cmath>
#include <assert.h>
#include <limits>

//#define USE_APPROX_COUNT
#define USE_EXACT_COUNT

typedef uint32_t contig_num_t;
#define MAX_CONTIG_NUM UINT32_MAX

typedef contig_num_t vertex_num_t;
typedef contig_num_t comp_num_t;
typedef uint32_t  rmer_count_t;
typedef uint8_t  kmer_len_t;
typedef uint32_t contig_edge_weight_t;

typedef uint64_t edge_weight_t;
typedef uint32_t read_num_t;
typedef uint32_t read_length_t;
typedef uint32_t read_count_t;

typedef uint64_t node_id_t;

typedef uint32_t kmer_count_node_t;


#ifdef USE_APPROX_COUNT
typedef uint16_t  kmer_count_t;
#define EXACT_THRESH 500
#define KMER_MAX_COUNT UINT16_MAX
#endif

#ifdef USE_EXACT_COUNT
typedef uint32_t  kmer_count_t;
#define KMER_MAX_COUNT UINT32_MAX
#endif

#define IMPOSSIBLE_READ_NUM (std::pow(2,sizeof(read_num_t)*8-1)-2)
#define IMPOSSIBLE_EDGE_WEIGHT_NUM (std::pow(2,sizeof(edge_weight_t)*8-1)-2)
#define IMPOSSIBLE_READ_LENGTH (std::pow(2,sizeof(read_length_t)*8)-1)

#define IMPOSSIBLE_NODE_ID (std::pow(2,sizeof(node_id_t)*8-3)-2)

#define NODE_ID_START 1000000000 //0
#define NODE_ID_END 1000000001  //1
#define NODE_ID_NORMAL_START 0 //2

#define IMPOSSIBLE_CONTIG_NUM (std::pow(2,sizeof(contig_num_t)*8-1)-2)
#define IMPOSSIBLE_COMP_NUM (std::pow(2,sizeof(comp_num_t)*8-2)-2)
#define IMPOSSIBLE_KMER_NUM (std::pow(2,sizeof(kmer_count_t)*8-2)-2)
#define IMPOSSIBLE_RMER_NUM (std::pow(2,sizeof(uint64_t)*8-3)-2)
/*
 * Info represents bit-field flag, six bits are used.
 * 1. The first bit describes if there is another kmer of suffix (k-1) common
 * characters relating to the prefix of this kmer. 1 for yes, 0 for no.
 * 2. The second bit describes if there is another kmer of prefix (k-1) common
 * characters relating to the suffix of this kmer. 1 for yes, 0 for no.
 * 3. The third, fourth bit indicates which "A,C,G,T" corresponding to the
 * prefix of this kmer relating to the first bit.
 * 4. The fifth, sixth bit indicates which "A,C,G,T" corresponding to the
 * suffix of this kmer relating to the second bit.
 */
struct Kmer_info {
    bool used;
    uint8_t info; // may refer to deleted kmer
    kmer_count_t count;
    contig_num_t contig;


    Kmer_info(): count(0), info(0), used(false), contig(IMPOSSIBLE_CONTIG_NUM){}
    Kmer_info(kmer_count_t cp): count(cp), info(0), used(false),
                contig(IMPOSSIBLE_CONTIG_NUM){}
    Kmer_info(kmer_count_t count_p, uint8_t info_p, bool up, contig_num_t contig_p
           ): count(count_p), info(info_p), used(up), contig(contig_p){}
};


#endif	/* SHC_TYPE_H */
