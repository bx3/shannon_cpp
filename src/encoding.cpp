#include "encoding.h"

uint8_t char_to_num [256] = {
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //   0
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //  20
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //  40
    255, 255, 255, 255, 255,   0, 255,   1, 255, 255, 255,   2, 255, 255, 255, 255, 255, 255, 255, 255, //  60
    255, 255, 255, 255,   3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,   0, 255,   1, //  80
    255, 255, 255,   2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,   3, 255, 255, 255, // 100
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 120
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 140
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 160
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 180
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 200
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 220
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255                      // 240
};

uint8_t num_to_char[4] = {'A','C','G','T'};
char xmer_buffer[32];

/**
 * As a general rule, say a list is encoded ATCG, where A is the first element
 * in an array and G is the last element, the corresponding bit format is 
 * [00....00(unused bit) 00(A) 11(T) 01(C) 10(G)]
 */

//length not include '\0'
void encode_kmer(const char * base, uint64_t *byte, uint8_t length)
{
    (*byte) = 0;
    for(int i=0; i<length; i++)
    {        
        *byte = (*byte)<< 2 | char_to_num[base[i]];
    }    
}

void decode_kmer(char *base, const uint64_t *byte, uint8_t length)
{
    uint64_t byte_copy = *byte;
    for (int i = length-1; i>=0; i--)
    {                
        base[i] = num_to_char[(uint8_t)B12 & byte_copy];    
        (byte_copy) >>=2;
    }
}

// for base less than or equal to 4, length exclude \0
uint8_t encode(const char * base, int length)
{
    uint8_t char_to_bin = 0x00;
    for(int i=0; i<length; i++)
    {        
        char_to_bin = char_to_bin<< 2 | char_to_num[base[i]];
    }
    return char_to_bin;
}

// for base less than or equal to 4, length exclude \0
uint8_t encode_in_place_helper(char * base, int length)
{
    uint8_t char_to_bin = 0x00;
    for(int i=0; i<length; i++)
    {        
        char_to_bin = char_to_bin<< 2 | char_to_num[base[i]];
    }
    return char_to_bin;
}

void decode(uint8_t num, int length, char *base)
{             
    for (int i=length-1; i>=0; i--)
    {                
        base[i] = num_to_char[(uint8_t)B12 & num];    
        num >>=2;
    }           
}

size_t encode_base_string(const char * base_string, uint8_t * byte_list, size_t base_length)
{    
    uint16_t byte_length = (uint16_t)ceil((base_length*1.0)/FOUR_BASE);
    
    int reminder_base_num = base_length % FOUR_BASE ; 
    
    int full_byte_num = 0;
    if (reminder_base_num > 0)
        full_byte_num = byte_length-1;
    else
        full_byte_num = byte_length;
    
    //printf("reminder %d, full_byte_num: %d\n",reminder_base_num,full_byte_num );
    
    int j = 0;
    
    for(int i = 0; i<full_byte_num; i++)
    {
        byte_list[j++] = encode(base_string+i*FOUR_BASE, FOUR_BASE);                 
        //printf("%d: %u\n", i, byte_list[j-1]);
        //print_four_base(base_string+i, FOUR_BASE);
        //printf("%u \t %d\n", byte_list[i], i);        
        //print_byte_list(byte_list, byte_length);
    }
    if (reminder_base_num > 0)
        byte_list[byte_length-1] = encode(base_string+FOUR_BASE*(byte_length-1), base_length%FOUR_BASE);
    return byte_length;
}

size_t encode_in_place(char * base_string, uint8_t * byte_list, size_t base_length)
{
    size_t byte_length = (size_t)ceil(static_cast<double>(base_length)/FOUR_BASE);
    
    int reminder_base_num = base_length % FOUR_BASE ; 
    
    int full_byte_num = 0;
    if (reminder_base_num > 0)
        full_byte_num = byte_length-1;
    else
        full_byte_num = byte_length;
    
    //printf("reminder %d, full_byte_num: %d\n",reminder_base_num,full_byte_num );
    
    int j = 0;
    uint8_t temp = 0;
    for(int i = 0; i<full_byte_num; i++)
    {
        //print_four_base(base_string+i*FOUR_BASE, FOUR_BASE);
        temp = encode_in_place_helper(base_string+i*FOUR_BASE, FOUR_BASE);        
        //printf("asa   %u\n", temp);
        byte_list[j++] = temp;
    }
    if (reminder_base_num > 0)
        byte_list[byte_length-1] = encode_in_place_helper(base_string+FOUR_BASE*(byte_length-1), base_length%FOUR_BASE);
    return byte_length;
}

void decode_byte_list(uint8_t * byte_list, char * base_string, size_t base_length)
{   
    uint8_t byte_length = (uint16_t)ceil(base_length*1.0/FOUR_BASE);
    
    int reminder_base_num = base_length % FOUR_BASE ; 
    //printf("byte_length is %u, %u \n", byte_length, reminder_base_num);
    int full_byte_num = 0;
    if (reminder_base_num > 0)
        full_byte_num = byte_length-1;
    else
        full_byte_num = byte_length;
    
    int j = 0;
    for(int i = 0; i<full_byte_num; i++)
    {
        decode(byte_list[i], FOUR_BASE, base_string+j);
        j+=4;                
    }
   
    if (reminder_base_num > 0)
        decode(byte_list[byte_length-1], reminder_base_num, base_string+j);               
}


void read_input(char* input, char* file_name)
{
    int c;
    int i = 0;
    FILE *file;
    file = fopen(file_name, "r");
    if (file) {
        while ((c = getc(file)) != EOF)
            input[i++] = c;            
        fclose(file);
    }
}

/**     Given a numeric kmer, get its first base in character*/
char decode_prefix_base(const uint64_t *byte, uint8_t length)
{
    return num_to_char[((char)((*byte) >> ((length-1)* BIT_PER_BASE))) & ((char)B12)];    
}

/**     Given a numeric kmer, get its first base in num*/
uint8_t decode_prefix_num(const uint64_t *byte, uint8_t length)
{
    return ((uint8_t)((*byte) >> ((length-1)* BIT_PER_BASE))) & ((uint8_t)B12);    
}

/**     Given a numeric kmer, get its last base in num*/
uint8_t decode_suffix_num(const uint64_t *byte)
{
    return (*byte) & ((uint64_t)B12);
}

uint64_t prepend_byte(const uint64_t * byte, uint64_t prefix, uint8_t length)
{
    return (*byte >> 2) | (prefix << (length-1)* BIT_PER_BASE);    
}

uint64_t append_byte(const uint64_t * byte, uint64_t prefix, uint8_t length)
{        
    return (((uint64_t)B1<<((length*2)))-(uint64_t)B1) &  ((*byte<<2)|prefix);
}

//asssume xmer > 4
void get_xmer_at_index(uint8_t * contig_start, size_t contig_base_len ,size_t index, 
                          uint8_t xmer_len, char * xmer_array)
{
    size_t xmer_start_byte_index = index/BASE_PER_BYTE;
    size_t xmer_end_byte_index = (index+xmer_len-1)/BASE_PER_BYTE;
    size_t total_byte = contig_base_len/BASE_PER_BYTE+1;
    
    uint8_t base_length =0;
    if (xmer_end_byte_index<total_byte-1 || (contig_base_len%BASE_PER_BYTE==0))
        base_length = (xmer_end_byte_index-xmer_start_byte_index+1)*BASE_PER_BYTE;
    else
        base_length = contig_base_len%BASE_PER_BYTE+(xmer_end_byte_index-xmer_start_byte_index)*BASE_PER_BYTE;
    
    memset(xmer_buffer,0,32);
    //printf("base_length %u\n", base_length);
    decode_byte_list(contig_start+xmer_start_byte_index, xmer_buffer, base_length);
    
    size_t index_in_buffer = index%BASE_PER_BYTE; 
    memcpy(xmer_array, xmer_buffer+index_in_buffer,  xmer_len);        
}


/**     Given a numeric kmer, get its last base in character*/
char decode_suffix_base(const uint64_t *byte)
{
    return num_to_char[(*byte) & ((uint64_t)B12)];
}

/**     
 *      combine multiple byte (<=4) into a 64bit primitive data
 *      Also notice that this implementation is consistent with
 *      encode_kmer function when less than 64bases<=>8byte, 
 *      a encode_kmer function gives the same output as
 *      encode_base_string followed with combine_byte.  
 */
uint64_t combine_byte(const uint8_t * byte, uint8_t length)
{
    uint64_t temp = 0;
    for (int i=0; i<length; i++)
    {        
        temp |= (uint64_t)byte[i] << ((length-i-1)*BIT_PER_BYTE);                
    }
    return temp;
}

void reverse_bit(uint8_t * byte, int len)
{
    uint8_t new_byte = (*byte)&(uint8_t)B12;    
    new_byte <<= 2;
    for(int i=0; i<len-2; i++)
    {        
        *byte >>= 2;
        new_byte |= ((uint8_t)B12&(*byte));
    }
    *byte = new_byte;
}

