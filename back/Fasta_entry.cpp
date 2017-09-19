
#include "Fasta_entry.h"
//#include "string_util.hpp"

void tokenize(const std::string& str, std::vector<std::string>& tokens, 
        const std::string& delimiters) {

    /* ************************************************************************************/
    /* from: http://oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html **/
    /**************************************************************************************/

    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }

    return;
}


// constructor
Fasta_entry::Fasta_entry(std::string header, std::string sequence) {
  
  // piece apart the header
  
  if (header[0] == '>') {
	header.erase(0,1);
  }

  
  std::vector<std::string> toks;
  std::string acc;
  
  tokenize(header, toks, " \t");
  
  if (toks.size() > 0) {
	acc = toks[0];
  }
  
  this->_header = header;
  this->_accession = acc;
  this->_sequence = sequence;
 

}

std::string Fasta_entry::get_accession() {
  return(this->_accession);
}

std::string Fasta_entry::get_header() {
  return(this->_header);
}

std::string Fasta_entry::get_sequence() {
  return(this->_sequence);
}
