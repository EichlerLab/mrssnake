/*

This program was created at:  Fri Feb 26 14:47:00 2016
This program was created by:  Brad Nelson


Contact: bnelsj@gmail.com

Organization: University of Washington

The MIT License (MIT)

Copyright (c) <2016> <Brad Nelson>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


*/

#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <inttypes.h>
#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/sam.h"

#define HTS_FMT_BAI 1

const char* const BAM_DNA_LOOKUP = "=ACMGRSVTWYHKDBN";
#define bam_seqi(s, i) ((s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf)

struct interval
{
  uint64_t start, end;
  bool operator < (const interval& other) const
  {
    return (start < other.start);
  }
};

struct options{
   std::string file;
   int part;
   int nparts;
   int chunk_size;
}globalOpts;

static const char *optString = "b:p:n:c";

//-------------------------------   OPTIONS   --------------------------------
int parseOpts(int argc, char** argv)
{
  int opt = 0;
  globalOpts.file = "NA";
  globalOpts.chunk_size = 36;
  opt = getopt(argc, argv, optString);
  while(opt != -1){
    switch(opt){
    case 'h':
      {
	std::cerr << "Useage" << std::endl;
      }
    case 'b':
      {
	globalOpts.file = optarg;
	break;
      }
    case 'p':
      {
	globalOpts.part = atoi(optarg);
	break;
      }
    case 'n':
      {
	globalOpts.nparts = atoi(optarg);
	break;
      }
    case 'c':
      {
	globalOpts.chunk_size = atoi(optarg);
	break;
      }
    case '?':
      {
	break;
      }
    }
    opt = getopt( argc, argv, optString ); 
  }
  return 1;
}

inline void initKstring(kstring_t * k){
  k->m = 0;
  k->l = 0;
  k->s = 0;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : nchunks in bam, partition number, total number of partitions

 Function does   : calculates start and end chunk to parse

 Function returns: half-open chunk interval for given partition

*/

std::vector<struct interval> get_chunk_range(std::vector<struct interval> chunks, int part, int nparts) {
	struct interval chunk_range;
	int range_size;
	int nchunks = chunks.size();

	range_size = chunks.size() / nparts;
	chunk_range.start = range_size * part;
	part == nparts-1 ? chunk_range.end = nchunks : chunk_range.end = range_size * part + range_size;
	
	//chunk_range.end = range_size * part + range_size;
	
	std::vector<struct interval> chunks_to_read(&chunks[chunk_range.start], &chunks[chunk_range.end]);

	return chunks_to_read;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : Vector of virtual offset intervals, bamfile

 Function does   : Iterate over vector, seek to start of interval and read to end

 Function returns: 0 if successful, 1 otherwise

*/
int get_reads(std::vector<interval> & chunks, 
	      std::string bamfile){
  
  BGZF* bam;
  bam = bgzf_open(bamfile.c_str(), "r");
  
  int index = 0;
  for(std::vector<interval>::iterator it = chunks.begin(); 
      it != chunks.end(); it++){
    
    index+= 1;
    
    if(index > 20){
      break;
    }
    
    if(bgzf_seek(bam, it->start, 0) == -1){
      std::cerr << "Error could not seek. " << std::endl;
      exit(1);
    }

    kstring_t seq = {0,0,NULL};

    
    int32_t r,refID, pos, lseq, nextRefID, nextPos, tlen;
    uint32_t bin_mq_nl, flag_nc, flag;
    
    
    
    std::cerr << sizeof(char) << std::endl;

    bgzf_read(bam, &r,         sizeof(r));
    bgzf_read(bam, &refID,     sizeof(refID));
    bgzf_read(bam, &pos,       sizeof(pos));
    bgzf_read(bam, &bin_mq_nl, sizeof(bin_mq_nl));
    bgzf_read(bam, &flag_nc,   sizeof(flag_nc));

    uint32_t test, cigar_ops;
    test = flag_nc >> 16;
	cigar_ops = flag_nc & 0xFFFF;
    

    bgzf_read(bam, &lseq,      sizeof(lseq));
    bgzf_read(bam, &nextRefID, sizeof(nextRefID));
    bgzf_read(bam, &nextPos,   sizeof(nextPos));
    bgzf_read(bam, &tlen,      sizeof(tlen));
    char readName[lseq];
    bgzf_read(bam, readName,   sizeof(char)*lseq);
	uint32_t cigar[cigar_ops];
    bgzf_read(bam, cigar,      sizeof(uint32_t)*cigar_ops);
    uint8_t seqa[(lseq+1)/2];
    bgzf_read(bam, seqa,       sizeof(uint8_t)*(lseq+1)/2);
    
    char seqASCII[lseq];
	for (int i = 0; i < (lseq+1)/2; ++i) seqASCII[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqa, i)];

    std::cerr << "remaining: " << r     << std::endl;
    std::cerr << "refID:     " << refID << std::endl;
    std::cerr << "pos:       " << pos   << std::endl;
    std::cerr << "lseq:      " << lseq   << std::endl;
    std::cerr << "nextrefID: " << nextRefID << std::endl;
    std::cerr << "nextPos  : " << nextPos << std::endl;
    std::cerr << "Readname  : " << readName << std::endl;
	std::cerr << "Flag: " << test << ", cigar_ops: " << cigar_ops << std::endl;
    std::cerr << "Cigar: " << cigar << std::endl;
    std::cerr << "Seq: " << seqASCII << std::endl;
    
    

  }
   bgzf_close(bam);
 
}
  
//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : options

 Function does   : get number of blocks in file, calculate start and end block from part and nparts

 Function returns: start and end positions

*/

void read_index(std::vector<interval> & chunks, 
		std::string & fn) {
 
 const char * file_name = fn.c_str();
  std::string index_name = fn + ".bai";
  long chunk_counter=0;
  FILE *ptr_myfile;
  //std::vector <struct interval> chunks;
  
  ptr_myfile = fopen(index_name.c_str(),"rb");
  
  if(!ptr_myfile) {
    std::cerr << "Unable to open file." << std::endl;
    exit(1);
  }
  
  char magic[4];
  fread(&magic, sizeof(char), 4, ptr_myfile);	
  printf("magic:%s\n",magic);
  
  int32_t n_ref;
  fread(&n_ref,sizeof(int32_t),1,ptr_myfile);
  printf("n_ref:%d\n",n_ref);
  
  int32_t nn;
  for(nn = 0; nn < n_ref; ++nn){
    
    int32_t n_bin;
    fread(&n_bin,sizeof(int32_t),1,ptr_myfile);
    printf(" n_bin:%lli\n",n_bin);
    
    int32_t i;
    for (i = 0; i < n_bin; ++i) {
      uint32_t bin; 
      fread(&bin,sizeof(uint32_t),1,ptr_myfile);
      
      int32_t n_chunk;
      fread(&n_chunk,sizeof(int32_t),1,ptr_myfile);
      printf("  n_chunk:%lli\n",n_chunk);
      
      int32_t j;
      for (j = 0; j < n_chunk; ++j) {
	uint64_t chunk_beg, chunk_end;
	interval tmp;
	fread(&chunk_beg,sizeof(uint64_t),1,ptr_myfile);
	fread(&chunk_end,sizeof(uint64_t),1,ptr_myfile);

	printf ( " offset %llu\n ", chunk_beg ) ;

	tmp.start = chunk_beg;
	tmp.end = chunk_end;
	chunks.push_back(tmp);		
      }
    }
    int32_t n_intv;
    fread(&n_intv,sizeof(int32_t),1,ptr_myfile);
    
    int32_t k; uint64_t ioffset;
    for(k = 0; k < n_intv; ++k){
      fread(&ioffset,sizeof(uint64_t),1,ptr_myfile);
    } 
  }

  uint64_t unmapped;
  fread(&unmapped, sizeof(uint64_t),1, ptr_myfile);
  std::cout << "unmapped " << unmapped << std::endl;

 
  fclose(ptr_myfile); 
}



//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  :

 Function does   :

 Function returns:

*/
/*
std::string chunk_read(std::string &read, int chunk_size)
{
	std::stringstream stream;
	int n_to_do = read->length / chunk_size;

    for(i = 0; i < n_to_do; i++){
        stream << ">0" << std::endl;
        stream << read->seq[i*chunk_size : i*chunk_size + chunk_size] << std::endl;
    }

    return stream;
}
*/
//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{
     
  int parse = parseOpts(argc, argv);
  if(parse != 1){
    std::cerr << "FATAL: unable to parse command line correctly." << std::endl;
    exit(1);
  }

 std::vector<interval> idx;

 read_index(idx, globalOpts.file);

 get_reads(idx, globalOpts.file);
    
  return 0;
}
