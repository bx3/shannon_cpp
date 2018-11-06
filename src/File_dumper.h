#ifndef FILE_DUMPER_H
#define	FILE_DUMPER_H

#include <memory>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>

#include <sys/ioctl.h>
#include <fcntl.h>

#include "log.h"

#define INIT_SIZE (1024*1024)
#define THRESH 0.95

#define FILE_MODE (S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)

#define NUM_OPEN_FILE 1000 //0
#define SRC_SIZE (2147483648)
//2147483648 536870912 104857600 500M  2G 100M

//#define LOG_FILER_DUMPER

typedef std::pair<int, uint64_t> comp_size_t;
struct File_size_sorter {
    bool operator()(const comp_size_t& cs1, const comp_size_t& cs2) const {
    return( cs1.second > cs2.second);
    }
};

struct File_chain_info {
    File_chain_info() :  num_file(0) {}
    File_chain_info(std::string filename_, int num_file_) :
                filename(filename_), num_file(num_file_) {}
    std::string filename;
    int num_file;
};

struct Src_man {
    Src_man() {}

    void setup(uint64_t total_num_char_, uint64_t mmap_sz_, std::string src_file_)
    {
        src_file = src_file_;
        header_ptr = (char*)src;
        curr_num_char = 0;
        local_file_num_char = 0;
        total_num_char = total_num_char_;
        mmap_sz = mmap_sz_;
        page_size = sysconf(_SC_PAGE_SIZE);
#ifdef LOG_FILER_DUMPER
        std::cout << "page size is " << page_size << std::endl;
#endif
    }

    void reset(void * src_, uint64_t mmap_sz_, char * header_ptr_)
    {
        src = (char*)src_;
        mmap_sz = mmap_sz_;
        header_ptr = header_ptr_;
    }

    void * src;  // used for free mmap
    int fdin;
    uint64_t mmap_sz;
    int page_size;

    uint64_t local_file_num_char;

    uint64_t curr_num_char;
    uint64_t total_num_char;

    char * header_ptr; // for locating read
    char * seq_ptr;

    std::string src_file;
};

struct File_man {
    File_man() : dst(NULL), offset(0), file_sz(1024), fdout(0),
                       num_resize(0), is_disk(true) {}
    File_man(bool is_disk_): is_disk(is_disk_) {}
    File_man(bool is_disk_, uint64_t file_sz_) : dst(NULL), offset(0),
                       file_sz(file_sz_), fdout(0), is_disk(is_disk_),
                       num_resize(0) {}

    File_man(bool is_disk_, uint64_t file_sz_, std::string filename_, int num_file_) :
                      dst(NULL), offset(0), is_disk(is_disk_),
                      file_sz(file_sz_), fdout(0),
                      num_resize(0),
                      chain_info(filename_, num_file_) {}
    void reset(uint64_t file_sz_, void * dst_, int fdout_)
    {
        file_sz = file_sz_;
        offset = 0;
        dst = dst_;
        fdout = fdout_;
    }

    bool is_disk;

    void * dst;
    uint64_t offset;
    uint64_t file_sz;
    int fdout;
    int num_resize;
    File_chain_info chain_info;

    std::string seq;
};

struct Single_dumper {

    Single_dumper() {}
    Single_dumper(uint64_t init_sz_, double increase_threshold_) :
               init_sz(init_sz_), increase_threshold(increase_threshold_)
    {
        //std::cout << "increase_threshold " << increase_threshold << std::endl;
        //std::cout << "init_size " << init_sz_ << std::endl;
        shc_log_info(shc_logname, "increase_threshold %f\n", increase_threshold);
        shc_log_info(shc_logname, "init_sz %lld\n", init_sz);

        resize_scale.push_back(2);
        resize_scale.push_back(2);
        resize_scale.push_back(2);
        resize_scale.push_back(2);
    }

    void setup_all_mmap_dst_files(
         std::string & file_prefix, std::string & suffix, int num_component)
    {

        num_comp = num_component;
        for(int i=0; i<num_component; i++)
        {
            std::string file_path = file_prefix + std::to_string(i) + suffix;
            std::set<int>::iterator it = disk_comps.find(i);
            if(it != disk_comps.end())
            {
                File_man file_man(true, init_sz);
                int fdout;
                void * dst;
                std::string partial_file_path(file_path + "_0");
                if ((fdout = open (partial_file_path.c_str(), O_RDWR | O_CREAT | O_TRUNC, FILE_MODE )) < 0)//edited here
                {   printf ("can't create for writing");
                    exit(1);
                }
                file_man.fdout = fdout;

                if(ftruncate(fdout, init_sz) == -1) {printf("ftruncate fail\n"); exit(1);}
                if ((dst = mmap (0, init_sz, PROT_WRITE | PROT_READ,
                                        MAP_SHARED, fdout, 0)) == (void*) -1)
                {perror ("mmap error for output\n"); exit(1); }
                file_man.dst = dst;

                file_man.chain_info.filename = file_path;
                file_man.chain_info.num_file++;

                comp_fileman_map.insert(std::make_pair(i, file_man));
            }
            else
            {
                File_man file_man(false);
                file_man.file_sz = init_sz;
                file_man.seq.reserve(file_man.file_sz);
                file_man.chain_info.filename = file_path;
                comp_fileman_map.insert(std::make_pair(i, file_man));
            }
        }
    }

    // return true if still lines to read
    inline bool get_src_line(char *& header_ptr,
                             char *& seq_ptr, int & seq_len)
    {
        if(is_increase_remap_src())
            re_mmap_src_file();

        if(src_man.curr_num_char < src_man.total_num_char)
        {
            int num_char_read = 0;
            char *ptr = (char*) memchr((void*)src_man.header_ptr, '\n',
                                src_man.total_num_char-src_man.curr_num_char);
            if (ptr == NULL)
            {
               printf("not found\n");
               exit(1);
            }
            //assign
            header_ptr = src_man.header_ptr;
            seq_ptr = ptr+1;
            src_man.curr_num_char += seq_ptr - header_ptr; // header len + 1
            src_man.local_file_num_char += seq_ptr - header_ptr;

            ptr = (char*) memchr((void*)seq_ptr, '\n',
                        src_man.total_num_char-src_man.curr_num_char);
            seq_len = ptr - seq_ptr;
            src_man.header_ptr = ptr + 1;
            src_man.curr_num_char += seq_len+1; // + 1 for \n
            src_man.local_file_num_char += seq_len+1;
            //shc_log_info(shc_logname, "%u\t%u\n", src_man.curr_num_char, src_man.local_file_num_char);
            return true;
        }
        else
        {
            return false;
        }
    }

    void setup_init_mmap_src_file(std::string input_read_path)
    {
        if ((src_man.fdin = open (input_read_path.c_str(), O_RDONLY)) < 0)
        {printf("can't open for reading"); exit(1);}
        if (fstat (src_man.fdin, &statbuf) < 0) {printf ("fstat error"); exit(1); }
#ifdef LOG_FILER_DUMPER
        std::cout << "statbuf.st_size is " << statbuf.st_size << std::endl;
#endif
        size_t mmap_size = 0;
        if(SRC_SIZE < statbuf.st_size)
            mmap_size = SRC_SIZE;
        else
            mmap_size = statbuf.st_size;

        if ((src_man.src = mmap (0, mmap_size, PROT_READ, MAP_SHARED, src_man.fdin, 0)) == (void*) -1)
        {printf ("mmap error for input"); exit(1); }
        src_man.setup(statbuf.st_size, mmap_size, input_read_path);
    }

    inline bool is_increase_remap_src()
    {
        //shc_log_info(shc_logname, "comp %d, len %lld, offset %lld, file_sz %lld\n",
        //                    comp_i, len, offsets[comp_i], files_sz[comp_i]);
        return static_cast<uint64_t>(increase_threshold * (SRC_SIZE))
                          <  src_man.local_file_num_char;
    }


    void re_mmap_src_file()
    {
#ifdef LOG_FILER_DUMPER
        std::cout << "remap src file" << std::endl;
#endif
        shc_log_info(shc_logname, "remap src file\n");

        if(munmap(src_man.src, src_man.mmap_sz) == -1){printf("munmap fail\n"); exit(1);}

        uint64_t num_extra_byte = src_man.curr_num_char % src_man.page_size;
        uint64_t offset = src_man.curr_num_char - num_extra_byte;

        size_t mmap_size = 0;
        if(SRC_SIZE < src_man.total_num_char - src_man.curr_num_char+offset)
            mmap_size = SRC_SIZE;
        else
            mmap_size = src_man.total_num_char - src_man.curr_num_char+offset;
        void * src;
        if ((src = mmap (src_man.src, mmap_size, PROT_READ, MAP_SHARED, src_man.fdin, offset)) == (void*) -1)
        {printf ("mmap error for input"); exit(1); }

        char * header_ptr = (char*)src + num_extra_byte;

        src_man.local_file_num_char = num_extra_byte;
        src_man.reset(src, mmap_size, header_ptr);
        shc_log_info(shc_logname, "finish remap src file\n");
    }

    void dump_file_increment_size(int i)
    {
        std::map<int, File_man>::iterator it = comp_fileman_map.find(i);
        if(it != comp_fileman_map.end())
        {
            File_man & file_man = it->second;
            std::string file_name(file_man.chain_info.filename + "_" +
                                std::to_string(file_man.chain_info.num_file++));
            int fdout;
            void * dst;

            if ((fdout = open (file_name.c_str(), O_RDWR | O_CREAT | O_TRUNC, FILE_MODE )) < 0)//edited here
            {   printf ("can't create for writing");
                exit(1);
            }

            uint64_t file_sz = file_man.seq.size();
            if(file_sz > 0 )
            {
                if(ftruncate(fdout, file_sz) == -1) {printf("ftruncate fail\n"); exit(1);}
                if ((dst = mmap (0, file_sz, PROT_WRITE | PROT_READ,
                                        MAP_SHARED, fdout, 0)) == (void*) -1)
                {perror ("mmap error for output\n"); exit(1); }
                memcpy(dst, &file_man.seq.at(0), file_sz);
                munmap(dst, file_sz);
            }
            close(fdout);
            file_man.seq.clear();
            uint64_t new_file_sz = get_new_sz(file_sz, i, file_man.num_resize);
            file_man.seq.reserve(new_file_sz);
            file_man.file_sz = new_file_sz;
            file_man.offset = 0;
        }
        else
        {
            std::cerr << "does not contain component " << i << std::endl;
            exit(1);
        }
    }

    void append_file_increment_size(int i)
    {
        std::map<int, File_man>::iterator it = comp_fileman_map.find(i);
        if(it != comp_fileman_map.end())
        {
            File_man & file_man = it->second;
            int fdout = file_man.fdout;
            void * dst = file_man.dst;
            //uint64_t offset = offsets[i];
            uint64_t file_sz = file_man.file_sz;
            if(ftruncate(fdout, file_man.offset) == -1) {printf("ftruncate fail\n"); exit(1);}
            // discard old mmap
            if(munmap(dst, file_sz) == -1){printf("munmap fail\n"); exit(1);}
            uint64_t new_file_sz = get_new_sz(file_sz, i, file_man.num_resize);
            close(fdout);

            std::string new_file(file_man.chain_info.filename + "_" +
                                std::to_string(file_man.chain_info.num_file++));
            //std::cout << "add file " << new_file << std::endl;
            if ((fdout = open (new_file.c_str(), O_RDWR | O_CREAT | O_TRUNC, FILE_MODE )) < 0)//edited here
            {   printf ("can't create for writing");
                exit(1);
            }

            if(ftruncate(fdout, new_file_sz) == -1) {printf("ftruncate fail\n"); exit(1);}
            if ((dst = mmap (0, new_file_sz, PROT_WRITE | PROT_READ,
                                    MAP_SHARED, fdout, 0)) == (void*) -1)
            {perror ("mmap error for output\n"); exit(1); }

            file_man.reset(new_file_sz, dst, fdout);
        }
        else
        {
            std::cerr << "does not contain component " << i << std::endl;
            exit(1);
        }
    }

    inline void mmap_write(uint64_t len, char * start, File_man & file_man)
    {
        //shc_log_info(shc_logname, "comp %d, len %lld, offset %lld, file_sz %lld\n",
        //                    comp_i, len, offsets[comp_i], files_sz[comp_i]);
        //memcpy(temp, header_ptr, len);
        //temp[len] = '\0';
        //shc_log_info(shc_logname, "pair %s", temp);
        //shc_log_info(shc_logname, "before write\n");

        void * write_ptr = (void * )((char*)file_man.dst + file_man.offset);
        memcpy(write_ptr, (void*)start, len);
        //shc_log_info(shc_logname, "after write\n");
        file_man.offset += len;
    }

    inline void reverse_complement(char * start_ptr, int len, char * write_ptr)
    {
        int j =0;
        for(int i=len-1; i>= 0; i--)
        {
            if (*(start_ptr+i)=='A')
                *(write_ptr + j++) = 'T';
            else if (*(start_ptr+i)=='T')
                *(write_ptr + j++) = 'A';
            else if (*(start_ptr+i)=='C')
                *(write_ptr + j++) = 'G';
            else if (*(start_ptr+i)=='G')
                *(write_ptr + j++) = 'C';
            else
            {
                //printf("unknown nucleotide %c\n", *(start_ptr+i));
                *(write_ptr + j++) = *(start_ptr+i);
                //exit(1);
            }
        }
    }

    inline char complement(char base)
    {
        if (base=='A')
            return 'T';
        else if ((base)=='T')
            return 'A';
        else if ((base)=='C')
            return 'G';
        else if ((base)=='G')
            return 'C';
        else
        {
            //printf("error unknown nucleotide %c\n", (base));
            return base;
        }

    }

    // len does not include \n
    inline void reverse_complement_local(char * start_ptr, int len)
    {
        char * s_ptr = start_ptr;
        char * e_ptr = start_ptr + len;
        while((s_ptr != e_ptr))
        {
            if (s_ptr!=--e_ptr)
            {
                char rc_base = complement(*s_ptr);
                *s_ptr = complement(*e_ptr);
                *e_ptr = rc_base;
                s_ptr++;
            }
            else
                *s_ptr = complement(*s_ptr);
        }
    }

    inline void reverse_complement_mem(char * start_ptr, int len, std::string & seq)
    {
        //std::string d_seq;

        for(int i=len-1; i>=0; i--)
        {
            seq += complement(*(start_ptr+i));
            //d_seq += complement(*(start_ptr+i));
        }
        //shc_log_info(shc_logname, "R %s\n", d_seq.c_str());
    }



    //len cover next line
    inline void mmap_reverse_write(uint64_t len, char * start, File_man & file_man)
    {
        char * write_ptr = ((char*)file_man.dst+file_man.offset);
        reverse_complement(start, len-1, write_ptr);
        //std::reverse_copy(start, start+len-1, write_ptr); // do not copy next line char

        *(write_ptr+len-1) = *(start+len-1);
        //memcpy(temp, start, len);
        //temp[len] = '\0';
        //shc_log_info(shc_logname, "comp %d reverse write %s\n", comp_i, temp);
        //shc_log_info(shc_logname, "after write\n");
        file_man.offset += len;
    }

    void deallocate_mmap()
    {
        close(src_man.fdin);
        munmap(src_man.src, src_man.mmap_sz);
        for(std::map<int, File_man>::iterator it=comp_fileman_map.begin();
                                    it!=comp_fileman_map.end(); it++)
        {
            int comp_i = it->first;
            File_man & file_man = it->second;
            if(file_man.is_disk)
            {
                if(ftruncate(file_man.fdout, file_man.offset) == -1) {printf("ftruncate fail\n"); exit(1);}
                close(file_man.fdout);
                munmap(file_man.dst, file_man.file_sz);
            }
        }
    }

    void deallocate_mem()
    {
        for(std::map<int, File_man>::iterator it=comp_fileman_map.begin();
                                    it!=comp_fileman_map.end(); it++)
        {
            int comp_i = it->first;
            File_man & file_man = it->second;
            if(!file_man.is_disk)
            {
                file_man.seq.clear();
            }
        }
    }

    inline uint64_t get_new_sz(uint64_t old_sz, int comp_i, int & resize_times)
    {
        resize_times++;
        if(resize_times < resize_scale.size())
        {
            return old_sz * resize_scale[resize_times];
        }
        else
        {
            return old_sz * resize_scale.back();
        }
    }

    inline bool is_increase_file_sz(uint64_t len, File_man & file_man)
    {
        //shc_log_info(shc_logname, "comp %d, len %lld, offset %lld, file_sz %lld\n",
        //                    comp_i, len, offsets[comp_i], files_sz[comp_i]);
        return static_cast<uint64_t>(increase_threshold *
                                     static_cast<double>(file_man.file_sz))
                          <  file_man.offset+len;
    }

    void set_disk_mem_components(std::set<int> & disk_comps_, std::set<int> & mem_comps_)
    {
        disk_comps = disk_comps_;
        mem_comps = mem_comps_;
    }

    inline void write_fasta(int comp_i, char * seq_ptr, char * header_ptr,
                                            int seq_len)
    {
        std::map<int, File_man>::iterator it = comp_fileman_map.find(comp_i);
        //std::string seq(seq_ptr, seq_len);
        //shc_log_info(shc_logname, "WF comp %d, seq %s\n", comp_i, seq.c_str());
        uint64_t len = seq_ptr - header_ptr + seq_len + 1; // include two next line
        if(it != comp_fileman_map.end())
        {
            File_man & file_man = it->second;
            if(file_man.is_disk)
            {
                if(is_increase_file_sz(len, file_man))
                {
                    append_file_increment_size(comp_i);
                }
                mmap_write(len, header_ptr, file_man);
            }
            else
            {
                if(is_increase_file_sz(len, file_man))
                {
                    dump_file_increment_size(comp_i);
                }
                //std::string header(header_ptr, len);
                //shc_log_info(shc_logname, "F %s\n", header.c_str());

                file_man.seq.append(header_ptr, len);
                file_man.offset = file_man.seq.size();
            }
        }
        else
        {
            std::cerr << "comp file manager does not exist" << std::endl;
            exit(1);
        }
    }

    inline void write_reverse_fasta(int comp_i, char * seq_ptr, char * header_ptr,
                                            int seq_len)
    {
        std::map<int, File_man>::iterator it = comp_fileman_map.find(comp_i);
        uint64_t len = seq_ptr - header_ptr + seq_len + 1;
        if(it != comp_fileman_map.end())
        {
            //std::string seq(seq_ptr, seq_len);
            //shc_log_info(shc_logname, "WR comp %d, seq %s\n", comp_i, seq.c_str());
            File_man & file_man = it->second;
            if(file_man.is_disk)
            {
                if(is_increase_file_sz(len, file_man))
                {
                    append_file_increment_size(comp_i);
                }
                mmap_write(seq_ptr-header_ptr, header_ptr, file_man); //write header
                mmap_reverse_write(seq_len+1, seq_ptr, file_man); //write reverse seq
            }
            else
            {
                if(is_increase_file_sz(len, file_man))
                {
                    dump_file_increment_size(comp_i);
                }

                //std::string header(header_ptr, seq_ptr-header_ptr); // include next line
                //shc_log_info(shc_logname, "%d %d header %s\n", seq_len, len, header.c_str());
                file_man.seq.append(header_ptr, seq_ptr-header_ptr);// += header;

                //std::string seq(seq_ptr, seq_len);
                //reverse_complement_local(seq_ptr, seq_len);
                //std::reverse(seq.begin(), seq.end());
                //file_man.seq.append(seq_ptr, seq_len);
                reverse_complement_mem(seq_ptr, seq_len, file_man.seq);
                file_man.seq += '\n';
            }
        }
        else
        {
            std::cerr << "comp file manager does not exist" << std::endl;
            exit(1);
        }
    }

    void finalize_dump_files()
    {
        Block_timer l_timer;
        start_timer(&l_timer);
        deallocate_mmap();
        /*
        std::cout << "start finalize disk read files" << std::endl;
        shc_log_info(shc_logname, "start finalize disk read files\n");
        for(std::set<int>::iterator it=disk_comps.begin();
                                    it!=disk_comps.end(); it++)
        {
            int comp_i = *it;
            std::map<int, File_man>::iterator fileman_it = comp_fileman_map.find(*it);
            if(fileman_it != comp_fileman_map.end())
            {
                File_man & file_man = fileman_it->second;
                std::string cmd("cat ");
                File_chain_info & chain_info = file_man.chain_info;
                for(int i=0; i<chain_info.num_file; i++)
                {
                    std::string file_path(chain_info.filename + "_" + std::to_string(i));
                    cmd += file_path + " ";
                }
                cmd += (" > " + chain_info.filename);
                Block_timer cat_timer;
                //start_timer(&cat_timer);
                run_command(cmd, false);
                //stop_timer(&cat_timer);
                for(int i=0; i<chain_info.num_file; i++)
                {
                    std::string rm_cmd("rm ");
                    std::string file_path(chain_info.filename + "_" +std::to_string(i));
                    rm_cmd += file_path;
                    run_command(rm_cmd, false);
                }

            }
            else
            {
                std::cerr << "cannot find such comp in finalize step" << std::endl;
                exit(1);
            }
        }
        std::cout << "finish finalize disk read files" << std::endl;
        stop_timer(&l_timer);
        log_stop_timer(&l_timer);
        */
#ifdef LOG_FILER_DUMPER
        std::cout << "start finalize mem read files" << std::endl;
#endif
        shc_log_info(shc_logname, "start finalize mem read files\n");
        start_timer(&l_timer);
        for(std::set<int>::iterator it=mem_comps.begin();
                                    it!=mem_comps.end(); it++)
        {
            int comp_i = *it;
            //shc_log_info(shc_logname, "dump mem for comp %d\n", comp_i);
            std::map<int, File_man>::iterator fileman_it = comp_fileman_map.find(*it);
            if(fileman_it != comp_fileman_map.end())
            {
                File_man & file_man = fileman_it->second;
                File_chain_info & chain_info = file_man.chain_info;
                int fdout;
                void * dst;
                std::string file_name(file_man.chain_info.filename + "_" +
                                    std::to_string(file_man.chain_info.num_file++));

                if ((fdout = open (file_name.c_str(), O_RDWR | O_CREAT | O_TRUNC, FILE_MODE )) < 0)//edited here
                {   printf ("can't create for writing");
                    exit(1);
                }

                uint64_t file_sz = file_man.seq.size();
                if(file_sz > 0 )
                {
                    if(ftruncate(fdout, file_sz) == -1) {printf("ftruncate fail\n"); exit(1);}
                    if ((dst = mmap (0, file_sz, PROT_WRITE | PROT_READ,
                                            MAP_SHARED, fdout, 0)) == (void*) -1)
                    {perror ("mmap error for output\n"); exit(1); }
                    memcpy(dst, &file_man.seq.at(0), file_sz);
                    munmap(dst, file_sz);
                }
                close(fdout);
                file_man.seq.clear();
            }
            else
            {
                std::cerr << "cannot find such comp in finalize step" << std::endl;
                exit(1);
            }
        }
#ifdef LOG_FILER_DUMPER
        std::cout << "finish finalize mem read files" << std::endl;
#endif
        //stop_timer(&l_timer);
        log_stop_timer(&l_timer);
        deallocate_mem();
    }

    int num_comp;
    uint64_t init_sz;

    std::vector<int> resize_scale;
    double increase_threshold;

    struct stat statbuf;

    Src_man src_man;

    std::map<int, File_man> comp_fileman_map;

    std::set<int> disk_comps;
    std::set<int> mem_comps;

    //debug
    char temp[10000];
};

struct FASTA_dumper {

    FASTA_dumper() {}
    FASTA_dumper(uint64_t init_sz_, double increase_threshold_,
                 bool has_single_, bool has_pair_) :
                single_dumper(init_sz_, increase_threshold_),
                pair1_dumper(init_sz_, increase_threshold_),
                pair2_dumper(init_sz_, increase_threshold_),
                has_single(has_single_), has_pair(has_pair_) {}

    void setup_dump_files_helper(std::vector<uint64_t> & component_size,
                bool is_single, std::set<int> & disk_comps,
                std::set<int> & mem_comps)
    {
        disk_comps.clear();
        if(is_single)
        {
            int num_allowed_comp = NUM_OPEN_FILE;
            int num_disk_comps = 0;
            if(num_component > num_allowed_comp)
            {
#ifdef LOG_FILER_DUMPER
                std::cout << num_allowed_comp << " files are opened for mmap dump reads" << std::endl;
#endif
                num_disk_comps = num_allowed_comp;
                std::vector<comp_size_t> comp_size_pair_array;
                for(int i=0; i<component_size.size(); i++)
                {
                    comp_size_t comp_size(i, component_size[i]);
                    comp_size_pair_array.push_back(comp_size);
                }
                File_size_sorter sorter;
                std::sort(comp_size_pair_array.begin(), comp_size_pair_array.end(), sorter);

                for(int i=0; i<num_component; i++)
                {
                    if(i < num_allowed_comp)
                    {
                        disk_comps.insert(comp_size_pair_array[i].first);
                        shc_log_info(shc_logname, "single comp %d is in disk\n",
                                        comp_size_pair_array[i].first);
                    }
                    else
                    {
                        mem_comps.insert(comp_size_pair_array[i].first);
                    }
                }
            }
            else
            {
#ifdef LOG_FILER_DUMPER
                std::cout << "all files are opened for mmap dump reads" << std::endl;
#endif
                num_disk_comps = num_component;
                for(int i=0; i<num_component; i++)
                {
                    disk_comps.insert(i);
                }
            }
        }
        else
        {
            int num_allowed_comp = NUM_OPEN_FILE/2;
            int num_disk_comps = 0;
            if(num_component > num_allowed_comp)
            {
#ifdef LOG_FILER_DUMPER
                std::cout << num_allowed_comp*2 << " files are opened for mmap dump reads" << std::endl;
#endif
                num_disk_comps = num_allowed_comp;
                std::vector<comp_size_t> comp_size_pair_array;
                for(int i=0; i<component_size.size(); i++)
                {
                    comp_size_t comp_size(i, component_size[i]);
                    comp_size_pair_array.push_back(comp_size);
                }
                File_size_sorter sorter;
                std::sort(comp_size_pair_array.begin(), comp_size_pair_array.end(), sorter);

                for(int i=0; i<num_component; i++)
                {
                    if(i < num_allowed_comp)
                    {
                        disk_comps.insert(comp_size_pair_array[i].first);
                        shc_log_info(shc_logname, "pair comp %d is in disk\n",
                                        comp_size_pair_array[i].first);
                    }
                    else
                    {
                        mem_comps.insert(comp_size_pair_array[i].first);
                        shc_log_info(shc_logname, "pair comp %d is in mem\n",
                                        comp_size_pair_array[i].first);
                    }
                }
            }
            else
            {
#ifdef LOG_FILER_DUMPER
                std::cout << "all files are opened for mmap dump reads" << std::endl;
#endif
                shc_log_info(shc_logname, "all pair files in disk\n");
                num_disk_comps = num_component;
                for(int i=0; i<num_component; i++)
                {
                    disk_comps.insert(i);
                }
            }
        }
    }

    void setup_dump_files (std::vector<uint64_t> & component_size)
    {
        shc_log_info(shc_logname, "Start setup_dump_files\n");
        num_component = component_size.size();
        if(has_single)
        {
            setup_dump_files_helper(component_size, true,
                        single_disk_comps, single_mem_comps);

            single_dumper.set_disk_mem_components(single_disk_comps, single_mem_comps);
        }

        if(has_pair)
        {
            setup_dump_files_helper(component_size, false,
                        pair_disk_comps, pair_mem_comps);
            pair1_dumper.set_disk_mem_components(pair_disk_comps, pair_mem_comps);
            pair2_dumper.set_disk_mem_components(pair_disk_comps, pair_mem_comps);
        }
        shc_log_info(shc_logname, "finish setup_dump_files\n");
    }

    void finalize_dump_files()
    {
        Block_timer dump_timer;
        if(has_single)
        {
#ifdef LOG_FILER_DUMPER
            std::cout << "Start single read dump" << std::endl;
#endif
            shc_log_info(shc_logname, "Start single read dump\n");
            start_timer(&dump_timer);
            single_dumper.finalize_dump_files();
#ifdef LOG_FILER_DUMPER
            std::cout << "finish dump single file" << std::endl;
            stop_timer(&dump_timer);
#endif
            log_stop_timer(&dump_timer);
        }
        if(has_pair)
        {
#ifdef LOG_FILER_DUMPER
            std::cout << "Start pair1 read dump" << std::endl;
#endif
            shc_log_info(shc_logname, "Start pair1 read dump\n");
            start_timer(&dump_timer);
            pair1_dumper.finalize_dump_files();
#ifdef LOG_FILER_DUMPER
            std::cout << "finish pair1 file" << std::endl;
            stop_timer(&dump_timer);
#endif
            log_stop_timer(&dump_timer);

#ifdef LOG_FILER_DUMPER
            std::cout << "Start pair2 read dump" << std::endl;
#endif
            shc_log_info(shc_logname, "Start pair2 read dump\n");
            start_timer(&dump_timer);
            pair2_dumper.finalize_dump_files();
#ifdef LOG_FILER_DUMPER
            std::cout << "finish pair2 file" << std::endl;
            stop_timer(&dump_timer);
#endif
            log_stop_timer(&dump_timer);
        }
        shc_log_info(shc_logname, "Finish finalize_dump_files\n");
    }

    int num_component;
    bool has_single;
    bool has_pair;

    struct Single_dumper single_dumper;
    std::set<int> single_disk_comps;
    std::set<int> single_mem_comps;

    struct Single_dumper pair1_dumper;
    struct Single_dumper pair2_dumper;
    std::set<int> pair_disk_comps;
    std::set<int> pair_mem_comps;
};

struct Kmer_man {
    Kmer_man() {}
    Kmer_man(bool is_disk_, std::string filename_) :
            is_disk(is_disk_), filename(filename_)
    {
        if(is_disk)
        {
            file_writer.open(filename_);
        }
    }
    bool is_disk;
    std::string filename;
    std::ofstream file_writer;
    std::string mem_writer;
};

struct KMER_dumper {
    KMER_dumper() {}

    void setup_dump_files (std::vector<uint64_t> & component_size, int num_allowed_comp_)
    {
        shc_log_info(shc_logname, "setup dump kmer\n");
        num_component = component_size.size();
        int num_allowed_comp = num_allowed_comp_;
        int num_disk_comps = 0;
        if(num_component > num_allowed_comp)
        {
            num_disk_comps = num_allowed_comp;
            std::vector<comp_size_t> comp_size_pair_array;
            for(int i=0; i<component_size.size(); i++)
            {
                comp_size_t comp_size(i, component_size[i]);
                comp_size_pair_array.push_back(comp_size);
            }
            File_size_sorter sorter;
            std::sort(comp_size_pair_array.begin(), comp_size_pair_array.end(), sorter);

            for(int i=0; i<num_component; i++)
            {
                if(i < num_allowed_comp)
                {
                    disk_comps.insert(comp_size_pair_array[i].first);
                    shc_log_info(shc_logname, "single comp %d is in disk\n",
                                    comp_size_pair_array[i].first);
                }
                else
                {
                    mem_comps.insert(comp_size_pair_array[i].first);
                }
            }
        }
        else
        {
            shc_log_info(shc_logname, "all single files in disk\n");
            num_disk_comps = num_component;
            for(int i=0; i<num_component; i++)
            {
                disk_comps.insert(i);
            }
        }
    }

    void setup_disk_dumper(std::string out_kmer_dir, std::string kmer_prefix)
    {
        for(comp_num_t i=0; i<num_component; i++)
        {
            std::set<int>::iterator it = disk_comps.find(i);
            std::string file_path = out_kmer_dir +
                            (kmer_prefix + std::to_string(i));
            if( it != disk_comps.end())
            {
                comp_kmer_map.emplace(i, Kmer_man(true, file_path));
            }
            else
            {
                comp_kmer_map.emplace(i, Kmer_man(false, file_path));
            }
        }
    }

    void dump_kmer(comp_num_t comp_i, char * base, kmer_count_t count)
    {
        std::string tab("\t");
        std::map<int, Kmer_man>::iterator it = comp_kmer_map.find(comp_i);
        if(it != comp_kmer_map.end())
        {
            Kmer_man & kmer_man = it->second;
            if(kmer_man.is_disk)
                kmer_man.file_writer << base << "\t" << count << std::endl;
            else
            {
                kmer_man.mem_writer +=
                            (base + tab + std::to_string(count) + "\n");
            }
        }
        else
        {
            std::cerr << "cannot find kmer dump" << std::endl;
            exit(1);
        }
    }

    void dump_read_prob(comp_num_t comp_i, double prob)
    {
        std::map<int, Kmer_man>::iterator it = comp_kmer_map.find(comp_i);
        if(it != comp_kmer_map.end())
        {
            Kmer_man & kmer_man = it->second;
            if(kmer_man.is_disk)
                kmer_man.file_writer << prob << std::endl;
            else
            {
                kmer_man.mem_writer +=
                            (std::to_string(prob) + "\n");
            }
        }
        else
        {
            std::cerr << "cannot find kmer dump" << std::endl;
            exit(1);
        }
    }

    void dump_read_features(comp_num_t comp_i, kmer_count_t * kmer_feat_array, int num_kmer)
    {
        std::string tab("\t");
        std::map<int, Kmer_man>::iterator it = comp_kmer_map.find(comp_i);
        if(it != comp_kmer_map.end())
        {
            Kmer_man & kmer_man = it->second;
            if(kmer_man.is_disk)
            {
                kmer_man.file_writer << kmer_feat_array[0];
                for(int i=1; i<num_kmer; i++)
                {
                    kmer_man.file_writer << tab  << (kmer_feat_array[i]);
                }
                kmer_man.file_writer << std::endl;
            }
            else
            {
                kmer_man.mem_writer +=  std::to_string(kmer_feat_array[0]);
                for(int i=1; i<num_kmer; i++)
                {
                    kmer_man.mem_writer += tab + std::to_string(kmer_feat_array[i]);
                }
                kmer_man.mem_writer += "\n";
            }
        }
        else
        {
            std::cerr << "cannot find kmer dump" << std::endl;
            exit(1);
        }
    }

    void dump_paired_read_features(comp_num_t comp_i,
                kmer_count_t * kmer_feat_array_1, int num_kmer1,
                kmer_count_t * kmer_feat_array_2, int num_kmer2)
    {
        std::string tab("\t");
        std::map<int, Kmer_man>::iterator it = comp_kmer_map.find(comp_i);
        if(it != comp_kmer_map.end())
        {
            Kmer_man & kmer_man = it->second;
            if(kmer_man.is_disk)
            {
                kmer_man.file_writer << kmer_feat_array_1[0];
                for(int i=1; i<num_kmer1; i++)
                {
                    kmer_man.file_writer << tab  << (kmer_feat_array_1[i]);
                }
                for(int i=0; i<num_kmer2; i++)
                {
                    kmer_man.file_writer << tab  << (kmer_feat_array_2[i]);
                }
                kmer_man.file_writer << std::endl;
            }
            else
            {
                kmer_man.mem_writer +=  std::to_string(kmer_feat_array_1[0]);
                for(int i=1; i<num_kmer1; i++)
                {
                    kmer_man.mem_writer += tab + std::to_string(kmer_feat_array_1[i]);
                }
                for(int i=0; i<num_kmer2; i++)
                {
                    kmer_man.mem_writer += tab + std::to_string(kmer_feat_array_2[i]);
                }
                kmer_man.mem_writer += "\n";
            }
        }
        else
        {
            std::cerr << "cannot find kmer dump" << std::endl;
            exit(1);
        }
    }

    void dump_read_count(comp_num_t comp_i, uint32_t count)
    {
        std::string tab("\t");
        std::map<int, Kmer_man>::iterator it = comp_kmer_map.find(comp_i);
        if(it != comp_kmer_map.end())
        {
            Kmer_man & kmer_man = it->second;
            if(kmer_man.is_disk)
                kmer_man.file_writer << count << std::endl;
            else
            {
                kmer_man.mem_writer +=
                            (std::to_string(count) + "\n");
            }
        }
        else
        {
            std::cerr << "cannot find kmer dump" << std::endl;
            exit(1);
        }
    }

    void write_mem_to_disk()
    {
        deallocate_mmap();

        for(std::set<int>::iterator it=mem_comps.begin();
                                    it!=mem_comps.end(); it++)
        {
            int comp_i = *it;
            std::map<int, Kmer_man>::iterator kmerman_it = comp_kmer_map.find(*it);
            if(kmerman_it != comp_kmer_map.end())
            {
                Kmer_man & kmer_man = kmerman_it->second;
                int fdout;
                void * dst;
                if ((fdout = open (kmer_man.filename.c_str(), O_RDWR | O_CREAT | O_TRUNC, FILE_MODE )) < 0)//edited here
                {   printf ("can't create for writing");
                    exit(1);
                }

                uint64_t file_sz = kmer_man.mem_writer.size();
                if(file_sz > 0 )
                {
                    if(ftruncate(fdout, file_sz) == -1) {printf("ftruncate fail\n"); exit(1);}
                    if ((dst = mmap (0, file_sz, PROT_WRITE | PROT_READ,
                                            MAP_SHARED, fdout, 0)) == (void*) -1)
                    {perror ("mmap error for output\n"); exit(1); }
                    memcpy(dst, &kmer_man.mem_writer.at(0), file_sz);
                    munmap(dst, file_sz);
                }
                close(fdout);
            }
            else
            {
                std::cerr << "cannot find such comp in finalize step" << std::endl;
                exit(1);
            }
        }
        deallocate_mem();
    }

    void deallocate_mem()
    {
        for(std::map<int, Kmer_man>::iterator it=comp_kmer_map.begin();
                                    it!=comp_kmer_map.end(); it++)
        {
            int comp_i = it->first;
            Kmer_man & kmer_man = it->second;
            if(!kmer_man.is_disk)
            {
                kmer_man.mem_writer.clear();
            }
        }
    }

    void deallocate_mmap ()
    {
        for(std::map<int, Kmer_man>::iterator it=comp_kmer_map.begin();
                                    it!=comp_kmer_map.end(); it++)
        {
            int comp_i = it->first;
            Kmer_man & kmer_man = it->second;
            if(kmer_man.is_disk)
            {
                kmer_man.file_writer.close();
            }
        }
    }

    void reset()
    {
        num_component = 0;
        comp_kmer_map.clear();
        disk_comps.clear();
        mem_comps.clear();
    }

    int num_component;
    std::map<int, Kmer_man> comp_kmer_map;
    std::set<int> disk_comps;
    std::set<int> mem_comps;
};


#endif
