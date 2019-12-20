/*
 * File:   shc_setting.h
 * Author: bx
 *
 * Created on January 5, 2018, 1:20 PM
 */

#ifndef RUN_TASKS_H
#define	RUN_TASKS_H

#include "Kmer_handler.h"
#include <unistd.h>
#include "Contig_graph_handler.h"
#include "encoding.h"
#include "local_file_structure.h"
#include "log.h"
#include <stdlib.h>
#include "json_parser.h"
#include "Multi_graph_handler.h"
#include "Sparse_flow_handler.h"
#include <sys/resource.h>
#include <stdio.h>
#include "shannon_C_seq_helper.h"
#include "Collect_reads.h"

//#include "genome_ref_filter.h"

//#define USE_MLOCK

void test_encoding_decoding();
void test_Kmer(Shannon_C_setting * setting);
void test_Contig_graph(Shannon_C_setting * setting);
void test_load_contig_graph(Shannon_C_setting * setting);
void test_load_dump(Shannon_C_setting * setting);
void test_pair_read_contig_graph(Shannon_C_setting * setting);
void test_load_pair_contig_graph(Shannon_C_setting * setting);
void test_seq_graph(Shannon_C_setting * setting);
void test_multi_seq_graph(Shannon_C_setting * setting);
void test_multithread_sparse_flow(Shannon_C_setting * setting);
void test_sparse_flow(Shannon_C_setting * setting);
void test_specific(Shannon_C_setting * setting);
void test_sorted_kmer_contig(Shannon_C_setting * setting);
void test_all(Shannon_C_setting * setting);
void test_multi_seq_sf(Shannon_C_setting * setting);
void test_assigning_reads_and_kmers(Shannon_C_setting * setting);
void test_evaluation(Shannon_C_setting * setting);
void test_custom_seq_graph(Shannon_C_setting * setting);
void sort_kmer(Shannon_C_setting * setting);
void test_filter_output(Shannon_C_setting * setting);
void test_load_contig_patition_and_multi_seq_sf(Shannon_C_setting * setting);
void test_sorted_kmer_contig_and_multi_seq_sf(Shannon_C_setting * setting);

void run_custom(int argc, char** argv, Shannon_C_setting setting);
void run_entire(int argc, char** argv);





#endif	/* RUN_TASKS_H */
