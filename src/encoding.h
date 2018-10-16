/*
 * File:   encoding.h
 * Author: bx
 *
 * Created on August 15, 2017, 4:37 PM
 */

#ifndef ENCODING_H
#define	ENCODING_H

#include <iostream>
#include <fstream>

#include <string.h>

#include "stdio.h"
#include <iostream>
#include "memory.h"
#include "stdint.h"
#include "math.h"

#include <map>

#include "log.h"

#define FOUR_BASE 4
#define BIT_PER_BASE 2
#define BIT_PER_BYTE 8
#define BASE_PER_BYTE 4
#define SHC_B78 0xC0        //set 7th, 8th bit up
#define SHC_B12 0x03        //set 1st, 2nd bit up
#define SHC_B1 0x01         //set 1st bit up
#define SHC_B2 0x02         //set 2nd bit up
#define SHC_B34 0x0C        //set 3th, 4th bit up
#define SHC_B56 0x30        //set 5th, 6th bit up
#define SHC_B1234 0x0F
#define SHC_B5678 0xF0
#define IS_BIT_I_SET(var,pos) ((var) & (1<<(pos)))

extern uint8_t char_to_num [256];
extern uint8_t num_to_char[4];

//length not include '\0'
inline void encode_kmer(const char * base, uint64_t *byte, uint8_t length)
{
    (*byte) = 0;
    for(int i=0; i<length; i++)
    {
        *byte = (*byte)<< 2 | char_to_num[base[i]];
    }
}

inline bool encode_kmer_check_base(const char * base, uint64_t *byte, uint8_t length)
{
    (*byte) = 0;
    for(int i=0; i<length; i++)
    {
        if(char_to_num[base[i]] >= 4)
            return false;
        *byte = (*byte)<< 2 | char_to_num[base[i]];
    }
    return true;
}

inline void encode_reverse_kmer(const char * base, uint64_t *byte, uint8_t length)
{
    (*byte) = 0;
    for(int i=length-1; i>=0; i--)
    {
        *byte = (*byte)<< 2 | char_to_num[base[i]];
    }
}

inline bool encode_reverse_kmer_check_base(const char * base, uint64_t *byte, uint8_t length)
{
    (*byte) = 0;
    for(int i=length-1; i>=0; i--)
    {
        if(char_to_num[base[i]] >= 4)
            return false;
        *byte = (*byte)<< 2 | char_to_num[base[i]];
    }
    return true;
}

inline void decode_kmer(char *base, const uint64_t *byte, uint8_t length)
{
    uint64_t byte_copy = *byte;
    for (int i = length-1; i>=0; i--)
    {
        base[i] = num_to_char[(uint8_t)SHC_B12 & byte_copy];
        (byte_copy) >>=2;
    }
}

inline void decode_reverse_kmer(char *base, const uint64_t *byte, uint8_t length)
{
    uint64_t byte_copy = *byte;
    for (int i = 0; i<length; i++)
    {
        base[i] = num_to_char[(uint8_t)SHC_B12 & byte_copy];
        (byte_copy) >>=2;
    }
}

inline uint8_t encode(const char * base, int length)
{
    uint8_t char_to_bin = 0x00;
    for(int i=0; i<length; i++)
    {
        char_to_bin = char_to_bin<< 2 | char_to_num[base[i]];
    }
    return char_to_bin;
}

inline void decode(uint8_t num, int length, char *base)
{
    for (int i=length-1; i>=0; i--)
    {
        base[i] = num_to_char[(uint8_t)SHC_B12 & num];
        num >>=2;
    }
}

inline uint64_t prepend_byte(const uint64_t * byte, const uint64_t num, uint8_t length)
{
    return (*byte >> 2) | (num << (length-1)* BIT_PER_BASE);
}

inline uint64_t append_byte(const uint64_t * byte, const uint64_t num, uint8_t length)
{
    return (((uint64_t)SHC_B1<<((length*2)))-(uint64_t)SHC_B1) &  ((*byte<<2)|num);
}

inline void complement_num(uint64_t *byte, uint8_t kmer_length)
{
    *byte = (~(*byte)) & ((((uint64_t)(0x01))<<(kmer_length*BIT_PER_BASE)) -
                                                            ((uint64_t)(0x01)));
}

inline bool is_info_ith_bit_set(uint8_t info, uint8_t i)
{
    if (i==1)
    {
        return (info & (uint8_t)(SHC_B1))>0;
    }
    else
    {
        return (info & (((uint8_t)(SHC_B1))<<(i-1)))>0;
    }
}

/**     Given a numeric kmer, get its first base in character*/
inline char decode_prefix_base(const uint64_t *byte, uint8_t length)
{
    return num_to_char[((char)((*byte) >> ((length-1)* BIT_PER_BASE))) & ((char)SHC_B12)];
}

inline uint8_t decode_prefix_num(const uint64_t *byte, uint8_t length)
{
    return ((char)((*byte) >> ((length-1)* BIT_PER_BASE))) & ((char)SHC_B12);
}


/**     Given a numeric kmer, get its last base in character*/
inline char decode_suffix_base(const uint64_t *byte)
{
    return num_to_char[(*byte) & ((uint64_t)SHC_B12)];
}

inline uint8_t decode_suffix_num(const uint64_t *byte)
{
    return (*byte) & ((uint64_t)SHC_B12);
}

uint64_t combine_byte(const uint8_t * byte, uint8_t length);
void get_xmer_at_index(uint8_t * contig_start, size_t contig_base_len ,
                       size_t index,uint8_t xmer_len, char * xmer_array);

uint8_t encode_in_place_helper(char * base, int length);
void decode(uint8_t num, int length, char *base);

size_t encode_base_string(const char * base_string, uint8_t * byte_list, size_t base_length);
size_t encode_in_place(char * base_string, uint8_t * byte_list, size_t base_length);

void decode_byte_list(uint8_t * byte_list, char * base_string, size_t base_length);

void read_input(char* input, char* file_name);

void complement_num(uint64_t *byte, uint8_t kmer_length);


#endif	/* ENCODING_H */
