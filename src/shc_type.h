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

typedef uint32_t contig_num_t;
typedef contig_num_t vertex_num_t;
typedef contig_num_t comp_num_t;
typedef uint32_t  kmer_count_t;
typedef uint16_t  rmer_count_t;
typedef uint8_t  kmer_len_t;
typedef uint16_t contig_edge_weight_t;

#define IMPOSSIBLE_CONTIG_NUM (std::pow(2,sizeof(contig_num_t)*8-2)-1)
#define IMPOSSIBLE_COMP_NUM (std::pow(2,sizeof(comp_num_t)*8-2)-2)
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
    kmer_count_t count;    
    uint8_t info; 
    bool used;
    contig_num_t contig;
    Kmer_info(): count(0), info(0), used(false), contig(IMPOSSIBLE_CONTIG_NUM){}
    Kmer_info(kmer_count_t cp): count(cp), info(0), used(false), 
                contig(IMPOSSIBLE_CONTIG_NUM){}
    Kmer_info(kmer_count_t count_p, uint8_t info_p, bool up, contig_num_t contig_p 
           ): count(count_p), info(info_p), used(false), contig(contig_p){}
};


#endif	/* SHC_TYPE_H */

