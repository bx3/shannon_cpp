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

void encode_kmer(const char * base, uint64_t *byte, uint8_t length);
void decode_kmer(char * base, const uint64_t *byte, uint8_t length);

bool is_info_ith_bit_set(uint8_t info, uint8_t i);

uint64_t prepend_byte(const uint64_t * byte, uint64_t prefix, uint8_t length);
uint64_t append_byte(const uint64_t * byte, uint64_t prefix, uint8_t length);

char decode_prefix_base(const uint64_t *byte, uint8_t length);
char decode_suffix_base(const uint64_t *byte);
uint8_t decode_prefix_num(const uint64_t *byte, uint8_t length);
uint8_t decode_suffix_num(const uint64_t *byte);

uint64_t combine_byte(const uint8_t * byte, uint8_t length);
void get_xmer_at_index(uint8_t * contig_start, size_t contig_base_len ,
                       size_t index,uint8_t xmer_len, char * xmer_array);

uint8_t encode(const char * base, int length);
uint8_t encode_in_place_helper(char * base, int length);
void decode(uint8_t num, int length, char *base);

size_t encode_base_string(const char * base_string, uint8_t * byte_list, size_t base_length);
size_t encode_in_place(char * base_string, uint8_t * byte_list, size_t base_length);

void decode_byte_list(uint8_t * byte_list, char * base_string, size_t base_length);

void read_input(char* input, char* file_name);

extern uint8_t char_to_num [256];
extern uint8_t num_to_char[4];


#endif	/* ENCODING_H */


