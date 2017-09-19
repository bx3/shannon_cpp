/* 
 * File:   File_reader.h
 * Author: bx
 *
 * Created on September 17, 2017, 10:26 PM
 */

#ifndef FILE_READER_H
#define	FILE_READER_H

#include <iostream>
#include <fstream>
#include <string>
#include "Fasta_entry.h"
#include <map>

using namespace std;

class Fasta_reader {
    
public:
    Fasta_reader(string filename); // constructor
    Fasta_reader(string filename, long start_reading, long end_reading); // constructor
    
    bool hasNext();  // prime for next 
    Fasta_entry getNext(); // retrieves next fasta entry, or NULL if EOF.
    
    map<string,string> retrieve_all_seqs_hash();
    
    
    
private:
    std::ifstream _filereader;
    //bool _hasNext;
    std::string _lastline;
    size_t end_reading; // optional file position to stop reading.

    size_t file_byte_pos;


    void _init_reader();

};

#endif	/* FILE_READER_H */

