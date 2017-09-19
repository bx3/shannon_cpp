/* 
 * File:   Fasta_entry.h
 * Author: bx
 *
 * Created on September 17, 2017, 10:30 PM
 */

#ifndef FASTA_ENTRY_H
#define	FASTA_ENTRY_H

#include <string>

class Fasta_entry {

public:

    Fasta_entry() {};
    Fasta_entry(std::string header, std::string sequence);

    std::string get_accession();

    std::string get_header();

    std::string get_sequence();


private:
    std::string _accession;
    std::string _header;
    std::string _sequence;
  
};

#endif	/* FASTA_ENTRY_H */

