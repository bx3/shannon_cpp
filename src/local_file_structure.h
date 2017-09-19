/* 
 * File:   local_file_structure.h
 * Author: bx
 *
 * Created on September 18, 2017, 1:10 PM
 */

#ifndef LOCAL_FILE_STRUCTURE_H
#define	LOCAL_FILE_STRUCTURE_H
#include <iostream>
#include <syscall.h>
#include <unistd.h>
#include <boost/filesystem.hpp>
#include <string.h>

void add_directory(boost::filesystem::path & dir_path);
void remove_directory(boost::filesystem::path & dir_path);

#endif	/* LOCAL_FILE_STRUCTURE_H */

