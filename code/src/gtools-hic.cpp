#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <map>
#include <string>
#include <list>
#include <vector>
#include <iostream>
#include <algorithm>
#include <sstream>
#include "core.h"
#include "gtools-intervals.h"
using namespace std;



// TODO
// ---------------
// * filter: more refined classification
// * [LATER] fragment reg file: check if non-overlapping (don't care about sorted)
// * [LATER] filter at the fragment level: too big, too small (--frag-min-size, --frag-max-size) [this needs to be done before creating the index]


//---------------------------------------------------------------------------------//
// Global variables & constants                                                    //
//---------------------------------------------------------------------------------//

const string PROGRAM = "gtools-hic";
const long int BUFFER_SIZE = 10000;
gsl_rng *RANDOM_GENERATOR;
const char *BIN_BITS = "17,20,23,26";


const int n_pair_types = 13;
enum PairType { UNPAIRED, UNMAPPED, MULTIHIT, SINGLE, DS_NO_FRAG, DS_SAME_FRAG, DS_TOO_CLOSE, DS_ACCEPTED_INTER, DS_ACCEPTED_INTRA, DS_DUPLICATE_INTER, DS_DUPLICATE_INTRA, DS_TOO_FAR, UNCLASSIFIED };
const char *PairTypeStr[] = { "unpaired", "unmapped", "multihit", "single-sided", "ds-no-fragment", "ds-same-fragment", "ds-too-close", "ds-accepted-inter", "ds-accepted-intra", "ds-duplicate-inter", "ds-duplicate-intra", "ds-too-far", "unclassified" };

typedef map<long int,bool> MapOfRead1;
typedef map<long int,MapOfRead1*> MapOfRead2;
typedef map<string,MapOfRead2*> MapOfRead;
MapOfRead dupMask;



//---------------------------------------------------------------------------------//
// COMMAND-LINE OPTIONS                                                            //
//---------------------------------------------------------------------------------//

bool VERBOSE;
bool HELP;
char *ALIGN_WORK_DIR;
int ALIGN_MIN_LEN;
int ALIGN_LEN_DIFF;
int BOWTIE_THREADS;
char *BOWTIE_PATH;
char *BOWTIE_INDEX;
char *ALIGNED_REG_FILE;
char *ENZYME_REG_FILE;
char *STATS_FILE;
double MAPQ;
long int FILTER_MIN_DISTANCE;
long int FILTER_MAX_OFFSET;
bool FILTER_DUPLICATES;
long int BIN_SIZE;
long int MAX_DIST;
bool ROTATE45;
char *GENOME_REG_FILE;
char *REF_REG_FILE;
char *OUTPUT_PREFIX;
bool PRINT_MATRIX_FORMAT;
bool SPLIT_MATRIX;
double NORM_CONSTANT;
double MIN_SCORE;
bool COL_LABELS;
char SEPARATOR; 



//-------InitCmdLine-----------
//
CmdLineWithOperations *InitCmdLine(int argc, char *argv[], int *next_arg)
{
  // initialize
  CmdLineWithOperations *cmd_line = new CmdLineWithOperations(); 
  cmd_line->SetProgramName(PROGRAM,VERSION);

  // set operations
  cmd_line->AddOperation("align", "[OPTIONS] READ1-FASTQ READ2-FASTQ", \
  "Iteratively aligns HiC-seq read pairs to reference genome using bowtie2.", \
  "* Input: FASTQ files\n\
  * Output: aligned reads in SAM format (same order as in fastq files)"\
  );
  
  cmd_line->AddOperation("classify", "[OPTIONS] <ALIGNED-READS>", \
  "Classifies and computes various metrics for HiC-seq aligned read pairs.", \
  "* Input: aligned reads in SAM format (sorted by read-id, at most one alignment per read)\n\
  * Output: tab-separated table"\
  );

  cmd_line->AddOperation("filter", "[OPTIONS] <ALIGNED-READS>", \
  "Filters HiC-seq aligned read pairs for common experimental artifacts.", \
  "* Input: aligned reads in SAM format (sorted by read-id, at most one alignment per read)\n\
  * Output: filtered read pairs in REG format"\
  );

  cmd_line->AddOperation("bin", "[OPTIONS] <FILTERED-READ-PAIRS>", \
  "Bins filtered read pairs to genomic bins of desired resolution.", \
  "* Input: filtered read pairs in REG format\n\
  * Output: binned read pairs"\
  );

  cmd_line->AddOperation("matrix", "[OPTIONS] <FILTERED-READ-PAIRS>", \
  "Create Hi-C count matrix.", \
  "* Input: filtered read pairs in REG format\n\
  * Output: contact matrix"\
  );

  cmd_line->AddOperation("convert", "[OPTIONS] <CONTACT-MATRIX>", \
  "Convert contact matrix into WashU Epigenome Browser format.", \
  "* Input: locus-labelled contact matrix\n\
  * Output: WashU Epigenome Browser format"\
  );

  map<string,int> min_files;
  
  // print summary of operations
  if (argc<2) { 
    cmd_line->OperationSummary("OPERATION [OPTIONS] <REGION-SET>", \
    "Pipeline for HiC-seq data analysis. For detailed description and list of options choose an operation and use the --help option.", true); 
    delete cmd_line; 
    exit(1); 
  }

  // current operation
  string op = argv[1];
  if (op[0]=='-') op = op.substr(1);		// ensure compatibility with previous version
  cmd_line->SetCurrentOperation(op);
  
  // common options
  cmd_line->AddOption("--help", &HELP, false, "help");
  cmd_line->AddOption("-h", &HELP, false, "help");
  cmd_line->AddOption("-v", &VERBOSE, false, "verbose mode");

  // Main options
  if (op=="align") {
    min_files["align"] = 2;
    cmd_line->AddOption("--work-dir", &ALIGN_WORK_DIR, "", "working directory (required)");
    cmd_line->AddOption("--min-len", &ALIGN_MIN_LEN, 30, "minimum truncated read length");
    cmd_line->AddOption("--len-diff", &ALIGN_LEN_DIFF, 10, "read truncation step");
    cmd_line->AddOption("-p", &BOWTIE_THREADS, 1, "number of threads for bowtie2 run");
    cmd_line->AddOption("--bowtie-path", &BOWTIE_PATH, "bowtie2", "full bowtie2 path (version>=2.1.0)");
    cmd_line->AddOption("--bowtie-index", &BOWTIE_INDEX, "genome/bowtie2.index/genome", "full bowtie2 index prefix path");
  }
  else if (op=="classify") {
    min_files["classify"] = 0;
    cmd_line->AddOption("-E", &ENZYME_REG_FILE, "", "enzyme fragments (BED/GFF/SAM/REG)");
    cmd_line->AddOption("--mapq", &MAPQ, 30, "minimum mapping quality (MAPQ)");
    cmd_line->AddOption("--min-dist", &FILTER_MIN_DISTANCE, 500, "miminum allowed distance between 5's of reads in read pair");
    cmd_line->AddOption("--max-offset", &FILTER_MAX_OFFSET, 500, "maximum allowed offset of 5's of reads from fragment ends");
  }
  else if (op=="filter") {
    min_files["filter"] = 0;
    cmd_line->AddOption("-E", &ENZYME_REG_FILE, "", "enzyme fragments (BED/GFF/SAM/REG)");
    cmd_line->AddOption("--mapq", &MAPQ, 30, "minimum mapping quality (MAPQ)");
    cmd_line->AddOption("--min-dist", &FILTER_MIN_DISTANCE, 500, "miminum allowed distance between 5's of reads in read pair");
    cmd_line->AddOption("--max-offset", &FILTER_MAX_OFFSET, 500, "maximum allowed offset of 5's of reads from fragment ends");
    cmd_line->AddOption("--filter-dups", &FILTER_DUPLICATES, false, "filter duplicate read pairs as PCR artifacts");
    cmd_line->AddOption("--stats", &STATS_FILE, "", "output statistics file (default=stderr)");
  }
  else if (op=="bin") {
    min_files["bin"] = 0;
    cmd_line->AddOption("--bin-size", &BIN_SIZE, 1000000, "genomic bin size");
    cmd_line->AddOption("-g", &GENOME_REG_FILE, "", "genome region file (BED/REG)");
    cmd_line->AddOption("--split-matrix", &SPLIT_MATRIX, false, "print output as matrix");
    cmd_line->AddOption("--matrix", &PRINT_MATRIX_FORMAT, false, "print output as matrix (overrides --split-matrix)");
  }
  else if (op=="matrix") {
    min_files["matrix"] = 0;
    cmd_line->AddOption("--bin-size", &BIN_SIZE, 5000, "genomic bin size (in nucleotides)");
    cmd_line->AddOption("--max-dist", &MAX_DIST, 0, "maximum distance between bins (in nucleotides; default = no restriction)");
    cmd_line->AddOption("--rotate45", &ROTATE45, false, "rotate matrix by 45 degrees (applicable if --max-dist > 0)");
    cmd_line->AddOption("-R", &REF_REG_FILE, "", "reference region file (BED/REG)");
    cmd_line->AddOption("-p", &OUTPUT_PREFIX, "", "output file prefix");
  }
  else if (op=="convert") {
    min_files["convert"] = 0;
    cmd_line->AddOption("--col-labels", &COL_LABELS, false, "input matrix has column labels");
    cmd_line->AddOption("-t", &SEPARATOR, ' ', "matrix element separator");
    cmd_line->AddOption("-c", &NORM_CONSTANT, 1.0, "normalization constant");
    cmd_line->AddOption("-min", &MIN_SCORE, 0.0, "score cutoff (values below this are set to zero)");
    cmd_line->AddOption("-d", &MAX_DIST, 0, "maximum distance between interacting loci (default = no limit)");
  }
  else {
    cerr << "Unknown operation '" << op << "'!\n";
    delete cmd_line;
    exit(1);
  }

  // process command line
  *next_arg = cmd_line->Read(argv+1,argc-1) + 1;
  if (HELP||(argc-*next_arg<min_files[op])) { cmd_line->OperationUsage(); delete cmd_line; exit(1); }
  
  return cmd_line;
}




//---MakeDirectory------
//
bool MakeDirectory(char *dir_name)
{
  struct stat st;
  if (stat(dir_name,&st)!=0) {
    if (VERBOSE) fprintf(stderr, "* Making directory '%s'\n", dir_name);
    mkdir(dir_name,0777);
    return true;
  }
  return false;
}


//---GetSAMFlagValue------
//
double GetSAMFlagValue(char *flagList, const char *flag)
{
  char *p = strstr(flagList,flag);
  if (p==NULL) return 0.0;
  char q[100];
  p = p+strlen(flag);
  for (int k=0; true; k++) { 
    q[k] = p[k]; 
    if ((p[k]=='\t')||(p[k]==0)) { q[k] = 0; break; }
  }
  double value = (double)atof(q);
  if (value<0) { fprintf(stderr, "Error: negative alignment score not allowed!\n"); exit(1); }
  return value;
}



//---IsReadPair------
//
bool IsReadPair(char *label1, char *label2)
{
  if (strcmp(label1,label2)==0) return true;                              // if read labels are identical, return true
  size_t n1 = strlen(label1);
  size_t n2 = strlen(label2);
  if (n1!=n2) return false;                                               // if different length, return false
  size_t n = n1;
  if (n<2) return false;                                                  // if too short, return false
  for (size_t k=0; k<n-1; k++) if (label1[k]!=label2[k]) return false;    // if there is a mismatch in the first n-1 characters, return false
  if (label1[n-2]!='.') return false;                                     // (n-1)-th character needs to be a '.'
  if ((label1[n-1]=='1')&&(label2[n-1]=='2')) return true;
  if ((label1[n-1]=='2')&&(label2[n-1]=='1')) return true;
  return false;
}


//---CLASS RegionPair------
//
template <class GenomicRegionType> class RegionPair
{
 public:
  RegionPair(FileBuffer *buffer, GenomicRegionSetIndex *refIndex);
  ~RegionPair();

  GenomicRegionType *r[2];
  GenomicRegion *r_matched[2];
  long int start_offset[2], stop_offset[2], min_offset[2];
  int n_frags[2];
  PairType p_type; 
};

//---constructor-------
//
template <class GenomicRegionType> RegionPair<GenomicRegionType>::RegionPair(FileBuffer *buffer, GenomicRegionSetIndex *refIndex)
{
  // read #1
  r[0] = new GenomicRegionType(buffer);
  
  // check if second read-id is present
  bool found;
  if (buffer->Get()!=NULL) {  
    char *p = StrCopy(buffer->Get());
    char *pp = p;
    char *read_id = GetNextToken(&pp,'\t');
    found = IsReadPair(r[0]->LABEL,read_id);
    free(p);
  }
  else found = false; 
  if (found==false) {
    p_type = UNPAIRED;
    r[1] = NULL;
    return;
  }

  // read #2
  r[1] = new GenomicRegionType(buffer);
  
  // reads are considered unique if they pass the minimum mappint quality threshold
  bool read1_uniq = r[0]->MAPQ >= MAPQ;
  bool read2_uniq = r[1]->MAPQ >= MAPQ;
  
  // assign read pair to corresponding enzyme fragments
  for (int q=0; q<=1; q++) {
    r_matched[q] = NULL;
    n_frags[q] = 0;
    start_offset[q] = stop_offset[q] = min_offset[q] = -1;
    if (strcmp(r[q]->I[0]->CHROMOSOME,"*")==0) continue;   // if unmapped, go to next read
    if (r[q]->I.size()!=1) r[q]->PrintError("alignment region is not contiguous!\n"); 
    GenomicInterval read5p = r[q]->I.front();
    read5p.ModifyPos("5p",0);       // assigned fragment will be the one containing the 5p of the read; therefore it will be unique
    r_matched[q] = refIndex->GetMatch(&read5p);
    if (r_matched[q]!=NULL) {
      n_frags[q] = 1; 
      r_matched[q]->I.front()->GetOffsetFrom(&read5p,"5p",true,&start_offset[q],&stop_offset[q]); 
      min_offset[q] = min(abs(start_offset[q]),abs(stop_offset[q]));
    }
  }

  // classify read pairs
  if ((strcmp(r[0]->I[0]->CHROMOSOME,"*")==0)&&(strcmp(r[1]->I[0]->CHROMOSOME,"*")==0)) p_type = UNMAPPED;     // both reads unmapped
  else if ((strcmp(r[0]->I[0]->CHROMOSOME,"*")==0)||(strcmp(r[1]->I[0]->CHROMOSOME,"*")==0)) p_type = SINGLE;    // one unmapped read
  else if ((read1_uniq==false)||(read2_uniq==false)) p_type = MULTIHIT;    // both reads mapped, but at least one is not mapped uniquely
  else if ((strcmp(r[0]->I[0]->CHROMOSOME,"*")!=0)&&(strcmp(r[1]->I[0]->CHROMOSOME,"*")!=0)) {   // both reads uniquely mapped
    // classify according to fragment assignment
    if ((n_frags[0]==0)||(n_frags[1]==0)) p_type = DS_NO_FRAG;   // at least one of the reads has no assigned fragment
    else {
      long int read_dist = abs(r[0]->I.front()->GetCoordinate("5p")-r[1]->I.front()->GetCoordinate("5p"));
      if (r_matched[0]==r_matched[1]) p_type = DS_SAME_FRAG;
      else if ((r[0]->I.front()->STRAND!=r[1]->I.front()->STRAND)&&(read_dist<FILTER_MIN_DISTANCE)) p_type = DS_TOO_CLOSE;
      else if ((min_offset[0]<=FILTER_MAX_OFFSET)&&(min_offset[1]<=FILTER_MAX_OFFSET)) {
        p_type = strcmp(r[0]->I[0]->CHROMOSOME,r[1]->I[0]->CHROMOSOME)==0?DS_ACCEPTED_INTRA:DS_ACCEPTED_INTER;
      }
      else if ((min_offset[0]>FILTER_MAX_OFFSET)||(min_offset[1]>FILTER_MAX_OFFSET)) p_type = DS_TOO_FAR;
    }
  }
  else p_type = UNCLASSIFIED;
}

//---destructor-------
template <class GenomicRegionType> RegionPair<GenomicRegionType>::~RegionPair()
{
  if (r[0]) delete r[0];
  if (r[1]) delete r[1];
}

//---ENDCLASS RegionPair-------



//------IsDuplicate------
//
bool IsDuplicate(GenomicInterval *r1, GenomicInterval *r2)
{
  string chr_pair = (string)r1->CHROMOSOME + r1->STRAND + (string)r2->CHROMOSOME + r2->STRAND;
  MapOfRead::iterator i_chr = dupMask.find(chr_pair);
  if (i_chr==dupMask.end()) {
    MapOfRead2 *m2 = dupMask[chr_pair] = new MapOfRead2();
    MapOfRead1 *m1 = (*m2)[r2->START] = new MapOfRead1();
    (*m1)[r1->START] = true;
    return false;
  }
  MapOfRead2 *m2 = i_chr->second;
  MapOfRead2::iterator im2 = m2->find(r2->START);
  if (im2==m2->end()) {
    MapOfRead1 *m1 = (*m2)[r2->START] = new MapOfRead1();
    (*m1)[r1->START] = true;
    return false;
  }
  MapOfRead1 *m1 = im2->second;
  MapOfRead1::iterator im1 = m1->find(r1->START);
  if (im1==m1->end()) {
    (*m1)[r1->START] = true;
    return false;
  }
  return true;
}




//------RunFilter------
//
void RunFilter(char **args, int argn, char *ref_reg_file, char *stats_file)
{
  // process args
  char *query_reg_file = argn==0?NULL:args[0];

  // check
  if (strlen(ref_reg_file)==0) { fprintf(stderr, "Error: need to supply enzyme fragment region file!\n"); exit(1); }
  
  // open region sets
  FileBuffer *QueryRegBuffer = CreateFileBuffer(query_reg_file,BUFFER_SIZE);
  GenomicRegionSet RefRegSet(ref_reg_file,BUFFER_SIZE,VERBOSE,true,true);
  GenomicRegionSetIndex RefIndex(&RefRegSet,BIN_BITS);

  // process read pairs
  Progress PRG("Filtering aligned read pairs...",1);
  unsigned long int n_reads = 0;
  unsigned long int n_class[n_pair_types];
  for (int i=0; i<n_pair_types; i++) n_class[i] = 0;
  for (char *qstr=QueryRegBuffer->Next(); QueryRegBuffer->Get()!=NULL; ) {
    if (qstr[0]=='@') { QueryRegBuffer->Next(); continue; }    // skip header lines
    RegionPair<GenomicRegionSAM> p(QueryRegBuffer,&RefIndex);
    if ((p.p_type==DS_ACCEPTED_INTER)||(p.p_type==DS_ACCEPTED_INTRA)) {
      GenomicInterval read1_5p = p.r[0]->I.front();
      read1_5p.ModifyPos("5p",0);
      GenomicInterval read2_5p = p.r[1]->I.front();
      read2_5p.ModifyPos("5p",0);
      bool duplicate = IsDuplicate(&read1_5p,&read2_5p);
      if (duplicate==true) {
        if (p.p_type==DS_ACCEPTED_INTER) p.p_type = DS_DUPLICATE_INTER;
        else if (p.p_type==DS_ACCEPTED_INTRA) p.p_type = DS_DUPLICATE_INTRA;
      }
      if ((FILTER_DUPLICATES==false)||(duplicate==false)) { 
        printf("%s\t", p.r[0]->LABEL);
        read1_5p.PrintInterval();
        printf(" ");
        read2_5p.PrintInterval();
        printf("\n");
      }
    }
    n_reads++;
    n_class[p.p_type]++;
    PRG.Check();
  }
  PRG.Done();

  // report
  FILE *fstats = strlen(stats_file)==0?stderr:fopen(STATS_FILE,"w");
  fprintf(fstats, "read-pairs\t%lu\t100%%\n", n_reads);
  for (int i=0; i<n_pair_types; i++) {
    fprintf(fstats, "%s\t", PairTypeStr[i]);
    fprintf(fstats, "%lu\t%.2f%%\n", n_class[i], (double)100*n_class[i]/n_reads);
  }
  if (strlen(stats_file)!=0) fclose(fstats);
  
  // cleanup
  delete QueryRegBuffer;
  for (MapOfRead::iterator i_chr=dupMask.begin(); i_chr!=dupMask.end(); i_chr++) {
    MapOfRead2 *m2 = i_chr->second;
    for (MapOfRead2::iterator i_m2=m2->begin(); i_m2!=m2->end(); i_m2++) delete i_m2->second;
    delete m2;
  }
}




//------RunClassify------
//
void RunClassify(char **args, int argn, char *ref_reg_file)
{
  // process args
  char *query_reg_file = argn==0?NULL:args[0];

  // check
  if (strlen(ref_reg_file)==0) { fprintf(stderr, "Error: need to supply enzyme fragment region file!\n"); exit(1); }
  
  // open region sets
  FileBuffer *QueryRegBuffer = CreateFileBuffer(query_reg_file,BUFFER_SIZE);
  GenomicRegionSet RefRegSet(ref_reg_file,BUFFER_SIZE,VERBOSE,true,true);
  GenomicRegionSetIndex RefIndex(&RefRegSet,BIN_BITS);

  // process read pairs
  Progress PRG("Classifying aligned read pairs...",1);
  unsigned long int n_reads = 0;
  unsigned long int n_class[n_pair_types];
  for (int i=0; i<n_pair_types; i++) n_class[i] = 0;
  printf("READ-LABEL");
  printf("\t");
  printf("READ1-SEQ\tREAD1-LOCUS\tREAD1-SIZE\tFRAG1-LABEL\tFRAG1-START-OFFSET\tFRAG1-STOP-OFFSET\tFRAG1-MIN-OFFSET\tFRAG1-LOCUS\tFRAG1-SIZE");
  printf("\t");
  printf("READ2-SEQ\tREAD2-LOCUS\tREAD2-SIZE\tFRAG2-LABEL\tFRAG2-START-OFFSET\tFRAG2-STOP-OFFSET\tFRAG2-MIN-OFFSET\tFRAG2-LOCUS\tFRAG2-SIZE");
  printf("\t");
  printf("DIST-NFRAGS\tDIST-FRAG-CENTERS");
  printf("\t");
  printf("PAIR-CHROMOSOMES\tPAIR-CLASSIFICATION\tPAIR-DUPLICATE\tPAIR-STRAND");
  printf("\n");
  for (char *qstr=QueryRegBuffer->Next(); QueryRegBuffer->Get()!=NULL; ) {
    if (qstr[0]=='@') { QueryRegBuffer->Next(); continue; }    // skip header lines
    
    RegionPair<GenomicRegionSAM> p(QueryRegBuffer,&RefIndex);

    // check if read pair is a duplicate
    bool is_duplicate = false;
    if ((p.p_type!=UNPAIRED)&&(p.p_type!=UNMAPPED)&&(p.p_type!=SINGLE)&&(p.p_type!=MULTIHIT)) {
      GenomicInterval read1_5p = p.r[0]->I.front();
      read1_5p.ModifyPos("5p",0);
      GenomicInterval read2_5p = p.r[1]->I.front();
      read2_5p.ModifyPos("5p",0);
      is_duplicate = IsDuplicate(&read1_5p,&read2_5p);
    }

    // print read info 
    printf("%s", p.r[0]->LABEL);     																	// READ-LABEL
    char strand[2] = { '.', '.' };
    for (int q=0; q<=1; q++) {
      printf("\t%s\t", p.r[q]->SEQ);         													// READ-SEQ
      p.r[q]->I.front()->PrintInterval();															// READ-LOCUS
      printf("\t%ld", p.r[q]->I.front()->GetSize());									// READ-SIZE
      if (p.r_matched[q]!=NULL) {
        printf("\t%s\t%ld\t%ld\t%ld\t", p.r_matched[q]->LABEL, p.start_offset[q], p.stop_offset[q], p.min_offset[q]);		// FRAG-LABEL, START-OFFSET, STOP-OFFSET, MIN-OFFSET
        p.r_matched[q]->I.front()->PrintInterval();										// FRAG-LOCUS
        printf("\t%ld", p.r_matched[q]->I.front()->GetSize());				// FRAG-SIZE
        strand[q] = p.r[q]->I.front()->STRAND;
      }
      else printf("\t*******\t*******\t*******\t*******\t*******\t*******");
    }

    // distances
    string chromdist;
    if ((p.r_matched[0]==NULL)||(p.r_matched[1]==NULL)) { chromdist = "N/A"; printf("\t*******\t*******"); }
    else if (strcmp(p.r_matched[0]->I.front()->CHROMOSOME,p.r_matched[1]->I.front()->CHROMOSOME)!=0) { chromdist = "interchrom"; printf("\tInf\tInf"); }
    else {
      chromdist = "intrachrom";
      printf("\t%ld", abs(atol(p.r_matched[0]->LABEL)-atol(p.r_matched[1]->LABEL)));																// DIST-NFRAGS
      printf("\t%ld", abs(p.r_matched[0]->I.front()->CalcDistanceFrom(p.r_matched[1]->I.front(),"c","c"))); 				// DIST-FRAG-CENTERS
    }
      
    // print classification results
    printf("\t%s", chromdist.c_str());																// PAIR-CHROMOSOMES
    printf("\t%s", PairTypeStr[p.p_type]);														// PAIR-CLASSIFICATION
    printf("\t%s", is_duplicate?"duplicate-pair":"unique-pair");			// PAIR-DUPLICATE
    printf("\tstrand=%c%c", strand[0], strand[1]);										// PAIR-STRAND
    printf("\n");
    
    // count
    n_reads++;
    n_class[p.p_type]++;
    PRG.Check();
  }
  PRG.Done();

  // report
  fprintf(stderr, "reads = %lu\n", n_reads);
  for (int i=0; i<n_pair_types; i++) {
    fprintf(stderr, "%s =", PairTypeStr[i]);
    fprintf(stderr, " %lu (%.2f%%)", n_class[i], (double)100*n_class[i]/n_reads);
    fprintf(stderr, "\n");
  }
  
  // cleanup
  delete QueryRegBuffer;
  for (MapOfRead::iterator i_chr=dupMask.begin(); i_chr!=dupMask.end(); i_chr++) {
    MapOfRead2 *m2 = i_chr->second;
    for (MapOfRead2::iterator i_m2=m2->begin(); i_m2!=m2->end(); i_m2++) delete i_m2->second;
    delete m2;
  }
}




//------RunAlign------
//
void RunAlign(char **args, int argn, char *work_dir, int min_len, int len_diff, int n_threads, char *bowtie_path, char *bowtie_index)
{
  // process args
  char *fastqFile1 = args[0];
  char *fastqFile2 = args[1];

  // determine read length
  FileBuffer *buffer = CreateFileBuffer(fastqFile1,BUFFER_SIZE); 
  buffer->Next();  // skip one line
  int read_len = strlen(buffer->Next());
  if (VERBOSE) fprintf(stderr, "* read length is %dnt\n", read_len);
  delete buffer;
  
  // initialize
  if (strlen(work_dir)==0) { fprintf(stderr, "Error: working directory must be specified!\n"); exit(1); }
  if (MakeDirectory(work_dir)==false) { fprintf(stderr, "Error: directory '%s' already exists!\n", work_dir); exit(1); }
  string bowtie_params = "--very-sensitive --reorder --sam-nohead -M 1";
  string fastq_file[2];
  string sam_file[2];
  fastq_file[0] = fastqFile1;
  fastq_file[1] = fastqFile2;
  for (int len=read_len; len>=min_len; len-=len_diff) {
    if (VERBOSE) fprintf(stderr, "Aligning reads (length=%d)...\n", len);
    for (int k=0; k<=1; k++) {
      // align
      stringstream ss_sam;
      ss_sam << work_dir << "/read" << k+1 << ".len=" << len << ".sam";
      sam_file[k] = ss_sam.str();
      stringstream cmd;
      cmd << bowtie_path << " -p " << n_threads << " " << bowtie_params << " -x " << bowtie_index << " -U " << fastq_file[k] << " -S " << sam_file[k];
      cerr << cmd.str() << '\n';
      system(cmd.str().c_str());
      // create new fastq files
      if (len-len_diff<min_len) continue;
      stringstream ss_q;
      ss_q << work_dir << "/read" << k+1 << ".len=" << len-len_diff << ".fastq";
      fastq_file[k] = ss_q.str();
      FileBuffer *buffer = CreateFileBuffer(sam_file[k].c_str(),BUFFER_SIZE); 
      FILE *Fout = fopen(fastq_file[k].c_str(),"w"); 
      Progress PRG("Processing SAM file...",1);
      for (buffer->Next(); buffer->Get()!=NULL; ) {
        GenomicRegionSAM r(buffer);
        GenomicInterval *i = r.I.front();
        if (strcmp(i->CHROMOSOME,"*")==0) fprintf(Fout, "@%s\n%.*s\n+%s\n%.*s\n", r.LABEL, len-len_diff, r.SEQ, r.LABEL, len-len_diff, r.QUAL);
        PRG.Check();
      }
      delete buffer;
      fclose(Fout);
      PRG.Done();
    }
  }

  // merge output sam files and print
  int n = (read_len-min_len)/len_diff+1;
  FileBuffer *samBuffer[2][n];
  for (int k=0; k<=1; k++) {
    for (int len=read_len,i=0; len>=min_len; len-=len_diff,i++) {
      stringstream ss_sam;
      ss_sam << work_dir << "/read" << k+1 << ".len=" << len << ".sam";
      samBuffer[k][i] = CreateFileBuffer(ss_sam.str().c_str(),BUFFER_SIZE); 
      samBuffer[k][i]->Next();
    }
  }
  Progress PRG2("Merging SAM files...",1);
  bool done = false;
  while (done==false) {
    for (int k=0; k<2; k++) {
      for (int i=0; i<n; i++) {
        char *s = samBuffer[k][i]->Get();
        if (s==NULL) { done = true; break; }
        GenomicRegionSAM r(s);
        samBuffer[k][i]->Next(); 
        if ((i==n-1)||(strcmp(r.I.front()->CHROMOSOME,"*")!=0)) { r.Print(); break; }
      }
    }
    PRG2.Check();
  }
  for (int k=0; k<=1; k++) for (int len=read_len,i=0; len>=min_len; len-=len_diff,i++) delete samBuffer[k][i];
  PRG2.Done();
}










//--------ProcessGenomeReg-----------
//
StringLIntMap *ProcessGenomeReg(char *genome_reg_file, long int bin_size, long int *n_bins, vector<string> &chr_names, map<string,long int> &chr_nbins, map<string,long int> &chr_k, map<string,long int> &chr_ends) 
{
  if ((genome_reg_file!=NULL)&&(strlen(genome_reg_file)>0)) {
    StringLIntMap *chrom_pos = new StringLIntMap();
    GenomicRegionSet RegSet(genome_reg_file,10000,false,false,true);
    long int line = 1;
    long int bin = 1;
    long int k = 0;
    for (GenomicRegion *r=RegSet.Get(); r!=NULL; r=RegSet.Next(),line++) {
      if (r->I.size()!=1) { cerr << "label = " << r->LABEL << '\n'; r->PrintError("genome regions should be single-interval regions!\n"); }
      if (r->I.front()->START!=1) { cerr << "label = " << r->LABEL << '\n'; r->PrintError("chromosome start coordinate should be 1!\n"); }
      string chr = r->I.front()->CHROMOSOME;
      if (chrom_pos->find(chr)==chrom_pos->end()) { 
        (*chrom_pos)[chr] = bin; 
        long int n = (long int)ceil((float)r->I.front()->GetSize()/bin_size); 
        bin += n;
        chr_names.push_back(chr);
        chr_nbins[chr] = n;
        chr_k[chr] = k;
        chr_ends[chr] = r->I.front()->STOP;
        k++;
      }
      else if ((*chrom_pos)[chr]!=r->I.front()->STOP) { cerr << "Error: chromosome " << chr << " appears in multiple lines in genome file '" << genome_reg_file << "' line " << line << "!\n"; exit(1); }
    }
    *n_bins = bin-1;
    return chrom_pos;
  }
  else {
    cerr << "Error: genome region file is necessary for this operation!\n"; exit(1);
    return NULL;
  } 
}



//------RunBin------
//
void RunBin(char **args, int argn, long int bin_size, char *genome_reg_file, bool print_matrix, bool split_matrix)
{
  // process args
  char *filteredRegFile = argn==0?NULL:args[0];

  // open region set
  FileBuffer *filteredRegBuffer = CreateFileBuffer(filteredRegFile,BUFFER_SIZE);

  // read genome reg file
  long int n_bins = 0;
  vector<string> chr_names;
  map<string,long int> chr_nbins;
  map<string,long int> chr_ends;
  map<string,long int> chr_k;  
  if (strlen(genome_reg_file)==0) { fprintf(stderr, "Error: genome region file is required for this operation!\n"); exit(1); }
  StringLIntMap *chrom_pos = ProcessGenomeReg(genome_reg_file,bin_size,&n_bins,chr_names,chr_nbins,chr_k,chr_ends);
  //for (StringLIntMap::iterator s=chrom_pos->begin(); s!=chrom_pos->end(); s++) cerr << s->first << '\t' << s->second << '\n';
  print_matrix = print_matrix&&(n_bins>0);
  long int **X = NULL;
  long int ***Y = NULL;
  if (print_matrix==true) {
    X = new long int*[n_bins];
    for (long int i=0; i<n_bins; i++) { X[i] = new long int[n_bins]; for (long int j=0; j<n_bins; j++) X[i][j] = 0; } 
  }
  else if (split_matrix==true) {
    long int n_chrom = chr_names.size();
    Y = new long int**[n_chrom];
    for (long int c=0; c<n_chrom; c++) {
      long int n_chrom_bins = chr_nbins[chr_names[c]];
      Y[c] = new long int*[n_chrom_bins];
      for (long int i=0; i<n_chrom_bins; i++) { Y[c][i] = new long int[n_chrom_bins]; for (long int j=0; j<n_chrom_bins; j++) Y[c][i][j] = 0; } 
    }
  }
  
  // process read pairs
  Progress PRG("Binning filtered reads...",1);
  unsigned long int n_discarded = 0;
  filteredRegBuffer->Next();   // get first line
  while (filteredRegBuffer->Get()!=NULL) {
    GenomicRegion r(filteredRegBuffer);
    if (r.I.size()!=2) { r.PrintError("each line should be a pair of aligned reads intervals"); exit(1); }
    char *bin1_chrom = r.I[0]->CHROMOSOME;
    char *bin2_chrom = r.I[1]->CHROMOSOME;
    long int bin1_offset = (r.I[0]->START-1)/bin_size+1;
    long int bin2_offset = (r.I[1]->START-1)/bin_size+1;
    StringLIntMap::iterator bin1_chrom_it = chrom_pos->find(bin1_chrom);
    StringLIntMap::iterator bin2_chrom_it = chrom_pos->find(bin2_chrom);
    if ((bin1_chrom_it!=chrom_pos->end())&&(bin2_chrom_it!=chrom_pos->end())) {
      long int i = bin1_chrom_it->second+bin1_offset-1;
      long int j = bin2_chrom_it->second+bin2_offset-1;
      if ((i>n_bins)||(j>n_bins)) { fprintf(stderr, "BUG!\n"); exit(1); }
      if (print_matrix==true) {
        X[i-1][j-1]++;
        if (i!=j) X[j-1][i-1]++;
      }
      else if (split_matrix==true) {
        if (strcmp(bin1_chrom,bin2_chrom)==0) {
          long int **X = Y[chr_k[bin1_chrom]];
          X[bin1_offset-1][bin2_offset-1]++;
          if (bin1_offset!=bin2_offset) X[bin2_offset-1][bin1_offset-1]++;
        }
        else n_discarded++;
      }
      else {
        long int start1 = bin_size*(bin1_offset-1)+1;
        long int stop1 = min(bin_size*bin1_offset,chr_ends[bin1_chrom]);
        long int start2 = bin_size*(bin2_offset-1)+1;
        long int stop2 = min(bin_size*bin2_offset,chr_ends[bin2_chrom]);
        printf("%s:%ld-%ld\t%s:%ld-%ld\n", bin1_chrom, start1, stop1, bin2_chrom, start2, stop2);
      }
    }
    else n_discarded++; 
    PRG.Check();
  }
  PRG.Done();

  // diagnostics
  if (VERBOSE) fprintf(stderr, "* discarded read pairs = %lu\n", n_discarded); 
  
  // print matrix
  if (print_matrix==true) {
    vector<string>::iterator chr_names_it = chr_names.begin();
    for (long int i=0,z=1; i<n_bins; i++) { 
      // print label
      long int w = chr_ends[*chr_names_it];
      cout << *chr_names_it << ":" << z << '-' << min(z+bin_size-1,w) << '\t'; 
      z += bin_size; 
      if (z>w) { ++chr_names_it; z = 1; } 
      // print values
      for (long int j=0; j<n_bins; j++) { printf("%ld", X[i][j]); printf("%c", j==n_bins-1?'\n':' '); }
    } 
  }
  else if (split_matrix==true) {
    vector<string>::iterator chr_names_it = chr_names.begin();
    long int n_chrom_bins = chr_nbins[*chr_names_it];
    long int c = 0;
    for (long int i=0,z=1; i<n_chrom_bins; i++) { 
      // print label
      long int w = chr_ends[*chr_names_it];
      cout << *chr_names_it << ":" << z << '-' << min(z+bin_size-1,w) << '\t'; 
      z += bin_size; 
      for (long int j=0; j<n_chrom_bins; j++) { printf("%ld", Y[c][i][j]); printf("%c", j==n_chrom_bins-1?'\n':' '); }
      if (z>w) {
        if (++chr_names_it==chr_names.end()) break; 
        ++c; 
        n_chrom_bins = chr_nbins[*chr_names_it]; 
        i = -1;
        z = 1; 
      } 
    } 
  }

  // cleanup
  if (chrom_pos!=NULL) delete chrom_pos;
  if (X!=NULL) delete [] X;
  if (Y!=NULL) { for (unsigned long int c=0; c<chr_names.size(); c++) delete [] Y[c]; delete Y; }
} 





// ###################################
//  CLASS HiCRefRegions
// ###################################
class HiCRefRegions
{
 public:
  HiCRefRegions(char *ref_reg_file, long int bin_size, long int max_dist, bool rotate45);
  ~HiCRefRegions();

  // methods
  void UpdateMatrix(GenomicRegion *r);
  void WriteOutput(char *file_prefix=NULL);
  bool CalcMatrixPosition(GenomicInterval *qint, long int *matrix_id, long int *matrix_offset);

  // data
  long int bin_size;
  long int max_dist;
  unsigned long int max_dist_bin;
  bool rotate45;
  bool match_gaps;
  bool ignore_strand;
  unsigned long int n_discarded;
  GenomicRegionSet *RefRegSet;
  GenomicRegionSet *RefIntSet;
  GenomicRegionSetIndex *RefIndex;
  unsigned long int n_intervals; 
  unsigned long int n_matrices; 
  unsigned long int ***X;               // array of contact matrices (ref-region x bin x bin)
  string *matrix_labels;                // matrix labels for each region
  unsigned long int *matrix_start;      // id of first interval in RefIntSet for each matrix
  unsigned long int *x_bins;            // number of row bins for each matrix
  unsigned long int *y_bins;            // number of column bins for each matrix
  unsigned long int *matrix_ids;        // matrix id number for each interval
  unsigned long int *matrix_offsets;    // matrix offset for each interval
};


// ###################################
//  Constructor
// ###################################
HiCRefRegions::HiCRefRegions(char *ref_reg_file, long int bin_size, long int max_dist, bool rotate45)
{
  // initialize
  this->bin_size = bin_size;
  this->max_dist = max(0L,max_dist);
  this->rotate45 = rotate45;
  this->match_gaps = false;
  this->ignore_strand = true;
  this->n_discarded = 0;
  this->max_dist_bin = (unsigned long int)ceil((float)max_dist/bin_size);
  
  // process reference region file
  if (VERBOSE) fprintf(stderr, "Processing reference regions...\n"); 
  this->RefRegSet = new GenomicRegionSet(ref_reg_file,10000,false,true,true);
  this->n_matrices = RefRegSet->n_regions;
  this->matrix_labels = new string[n_matrices];
  this->matrix_start = new unsigned long[n_matrices];
  this->x_bins = new unsigned long int[n_matrices];
  this->y_bins = new unsigned long int[n_matrices];
  this->n_intervals = 0;
  for (unsigned long int k=0; k<n_matrices; k++) n_intervals += RefRegSet->R[k]->I.size();
  this->RefIntSet = new GenomicRegionSet(n_intervals,RefRegSet->format,false);
  if (VERBOSE) fprintf(stderr, "* Found %lu regions and %lu intervals.\n", n_matrices, n_intervals); 
  this->matrix_offsets = new unsigned long int[n_intervals];
  this->matrix_ids = new unsigned long int[n_intervals];
  map<string,bool> label_map;
  for (unsigned long int k=0,q=0; k<n_matrices; k++) {
    GenomicRegion *rreg = RefRegSet->R[k];
    x_bins[k] = 0;
    if (label_map.find(rreg->LABEL)!=label_map.end()) { fprintf(stderr, "Error: reference region labels should be unique!\n"); exit(1); }
    label_map[rreg->LABEL] = true;
    matrix_labels[k] = rreg->LABEL;
    matrix_start[k] = q;
    for (GenomicIntervalSet::iterator i=rreg->I.begin(); i!=rreg->I.end(); i++,q++) {
      RefIntSet->R[q] = new GenomicRegion(rreg->LABEL,(*i));
      RefIntSet->R[q]->n_line = q;
      matrix_ids[q] = k;
      matrix_offsets[q] = x_bins[k];
      x_bins[k] += ((*i)->GetSize()-1)/bin_size + 1;
    }
    if (max_dist<=0) y_bins[k] = x_bins[k];
    else if (rotate45==false) y_bins[k] = 2*max_dist_bin+1;
    else y_bins[k] = max_dist_bin+1;
    if (VERBOSE) fprintf(stderr, "* Region %s row/column bins = %lu/%lu\n", rreg->LABEL, x_bins[k], y_bins[k]); 
  }

  // create index for reference intervals
  this->RefIndex = new GenomicRegionSetIndex(RefIntSet,BIN_BITS);

  // initialize contact matrix
  if (VERBOSE) fprintf(stderr, "Initializing contact matrices...\n"); 
  X = new unsigned long int**[n_matrices];
  for (unsigned long int k=0; k<n_matrices; k++) {
    X[k] = new unsigned long int*[x_bins[k]];
    if (X[k]==NULL) { fprintf(stderr, "Error: out of memory!\n"); exit(1); }
    for (unsigned long int i=0; i<x_bins[k]; i++) { 
      X[k][i] = new unsigned long int[y_bins[k]]; 
      if (X[k][i]==NULL) { fprintf(stderr, "Error: out of memory!\n"); exit(1); }
      for (unsigned long int j=0; j<y_bins[k]; j++) X[k][i][j] = 0; 
    } 
  }
  if (VERBOSE) fprintf(stderr, "Initialization completed.\n"); 
}



// ###################################
//  Destructor
// ###################################
HiCRefRegions::~HiCRefRegions()
{
//  delete matrix_labels;
  delete matrix_start;
  delete RefRegSet;
  delete RefIntSet;
  delete RefIndex;
  delete x_bins;
  delete y_bins;
  delete matrix_ids;
  delete matrix_offsets;
  for (unsigned long int k=0; k<n_matrices; k++) delete [] X[k];
  delete X;
}



// ###################################
//  CalcMatrixPosition
// ###################################
bool HiCRefRegions::CalcMatrixPosition(GenomicInterval *qint, long int *matrix_id, long int *matrix_offset)
{
  GenomicRegion *ireg = RefIndex->GetOverlap(qint,match_gaps,ignore_strand);   // TODO: generalize to deal with multi-hits?
  if (ireg==NULL) {
    *matrix_id = *matrix_offset = -1;
    return(false);
  }
  unsigned long int q = ireg->n_line;
  *matrix_id = matrix_ids[q];
  *matrix_offset = matrix_offsets[q] + (qint->START - ireg->I[0]->START) / bin_size;
  return(true);
}



// ###################################
//  UpdateMatrix
// ###################################
void HiCRefRegions::UpdateMatrix(GenomicRegion *r)
{
  if (r->I.size()!=2) { r->PrintError("each line should be a pair of aligned reads intervals"); exit(1); }
  long int id1, id2;
  long int i, j;
  bool b1 = CalcMatrixPosition(r->I[0],&id1,&i);
  bool b2 = CalcMatrixPosition(r->I[1],&id2,&j);
  if (b1 && b2) {
    if (id1==id2) {
      unsigned long int k = (unsigned long int)id1;
      if (max_dist<=0) {   // no restriction
        if ((k>=n_matrices)||(i>=(long int)x_bins[k])||(j>=(long int)y_bins[k])) { fprintf(stderr, "Error: matrix index error, this should be a bug!\n"); exit(1); } 
        X[k][i][j]++;
        if (i!=j) X[k][j][i]++;
      }
      else if (rotate45==false) {
        unsigned long int d = abs(j-i);
        if (d<=max_dist_bin) { 
          if ((k>=n_matrices)||(max(i,j)>=(long int)x_bins[k])||(max_dist_bin+d>=y_bins[k])) { fprintf(stderr, "Error: matrix index error, this should be a bug!\n"); exit(1); } 
          X[k][i][max_dist_bin+j-i]++;
          if (d!=0) X[k][j][max_dist_bin+i-j]++;
        }
      }
      else {
        long int j2 = abs(i-j);
        long int i2 = (i+j)/2;
        if ((k>=n_matrices)||(i2>=(long int)x_bins[k])) { fprintf(stderr, "Error: matrix index error, this should be a bug!\n"); exit(1); } 
        if (j2<(long int)y_bins[k]) X[k][i2][j2]++;
      }
    }
    else n_discarded++; 
  }
  else n_discarded++; 
}


// ###################################
//  WriteOutput
// ###################################
void HiCRefRegions::WriteOutput(char *file_prefix)
{
  // save matrices
  for (unsigned long int k=0; k<n_matrices; k++) {
    // open output file
    string fname = file_prefix + matrix_labels[k] + ".tsv";
    string msg = "Storing " + fname + " matrix...";
    FILE *fout = fopen(fname.c_str(),"w"); 
    if (fout==NULL) { fprintf(stderr, "Error: cannot open file '%s' for writing!\n", fname.c_str()); exit(1); }

    // initialize
    unsigned long int z;        // RefIntSet index
    GenomicInterval *rint;      // RefIntSet interval
    long int start;             // interval start position

    // print column labels (for matrices that are not distance-restricted)
    if (max_dist<=0) {    
      rint = RefIntSet->R[z=matrix_start[k]]->I[0];
      start = rint->START;
      for (unsigned long int i=0; i<y_bins[k]; i++) { 
        if (start>rint->STOP) { rint = RefIntSet->R[++z]->I[0]; start = rint->START; } 
        fprintf(fout, "\t%s:%ld-%ld", rint->CHROMOSOME, start, min(start+bin_size-1,rint->STOP));
        start += bin_size;
      }
      fprintf(fout, "\n");
    }

    // print matrix
    Progress PRG(msg.c_str(),x_bins[k]);
    rint = RefIntSet->R[z=matrix_start[k]]->I[0];
    start = rint->START;
    for (unsigned long int i=0; i<x_bins[k]; i++) { 
      // print row label
      if (start>rint->STOP) { rint = RefIntSet->R[++z]->I[0]; start = rint->START; } 
      fprintf(fout, "%s:%ld-%ld\t", rint->CHROMOSOME, start, min(start+bin_size-1,rint->STOP));
      start += bin_size;
      // print values
      for (unsigned long int j=0; j<y_bins[k]; j++) { fprintf(fout, "%ld", X[k][i][j]); fprintf(fout, "%c", j==y_bins[k]-1?'\n':'\t'); }
      PRG.Check();
    }

    // done
    fclose(fout);
    PRG.Done();
  }
}


// ###################################
//  END CLASS HiCRefRegions
// ###################################




//------RunMatrix------
//
void RunMatrix(char **args, int argn, long int bin_size, long int max_dist, bool rotate45, char *ref_reg_file, char *output_prefix)
{
  // check parameters
  if ((ref_reg_file==NULL)||(strlen(ref_reg_file)==0)) {
    cerr << "Error: genome region file is necessary for this operation!\n"; 
    exit(1);
  }
  
  // initialize reads and reference regions
  char *filteredRegFile = argn==0?NULL:args[0];
  FileBuffer *filteredRegBuffer = CreateFileBuffer(filteredRegFile,BUFFER_SIZE);
  HiCRefRegions ref(ref_reg_file,bin_size,max_dist,rotate45);
  
  // overlaps & binning
  Progress PRG("Assigning reads pairs to bin pairs...",1);
  filteredRegBuffer->Next();                  // get first line
  while (filteredRegBuffer->Get()!=NULL) {
    GenomicRegion qreg(filteredRegBuffer);       // get filtered read pair
    ref.UpdateMatrix(&qreg);
    PRG.Check();
  }
  if (VERBOSE) fprintf(stderr, "* Discarded reads = %lu.\n", ref.n_discarded);
  PRG.Done();
  
  // write contact matrices to output
  ref.WriteOutput(output_prefix);
  
} 



//------GetLabelFields-------
//
string GetLabelFields(string *lab, string &coord)
{
  size_t p = lab->find(':');
  string chrom = lab->substr(0,p);
  coord = lab->substr(p+1);
  coord[coord.find('-')] = '\t';
  return chrom;
}


//------RunConvert------
//
void RunConvert(bool col_labels, char separator, double norm_constant, double min_score, long int max_dist)
{
  bool row_labels = true;
  MatrixTemplate<double> M(NULL,VERBOSE,separator,row_labels,col_labels);
  Progress PRG("Converting...",M.n_rows);
  unsigned long int id = 1;
  for (long int r=0; r<M.n_rows; r++) {
    string coord1;
    string chrom1 = GetLabelFields(&M.row_labels[r],coord1);
    for (long int c=0; c<M.n_cols; c++) {
      M.val[r][c] /= norm_constant; 
      if (fabs(M.val[r][c])>min_score) {
        string coord2;
        string chrom2 = GetLabelFields(&M.row_labels[c],coord2);
        char flag = chrom1!=chrom2?'.':(c>=r?'+':'-');
        long int dist = abs(atol(coord1.c_str())-atol(coord2.c_str()));
        if ((max_dist<=0)||(chrom1!=chrom2)||(dist<=max_dist)) cout << chrom1 << '\t' << coord1 << '\t' << M.row_labels[c] << ',' << M.val[r][c] << '\t' << id++ << '\t' << flag << '\n';
      }
    }
    PRG.Check();
  }
  PRG.Done();
}





//---------------------------------------------------------------------------------//
// MAIN	                                                                           //
//---------------------------------------------------------------------------------//
int main(int argc, char* argv[]) 
{
  // process command-line arguments
  int next_arg;
  CmdLineWithOperations *cmd_line = InitCmdLine(argc,argv,&next_arg); 
  _MESSAGES_ = VERBOSE;

  // initilize random generator
  RANDOM_GENERATOR = InitRandomGenerator(getpid()+time(NULL));

  // run 
  if (cmd_line->current_cmd_operation=="align") RunAlign(&argv[next_arg],argc-next_arg,ALIGN_WORK_DIR,ALIGN_MIN_LEN,ALIGN_LEN_DIFF,BOWTIE_THREADS,BOWTIE_PATH,BOWTIE_INDEX);
  else if (cmd_line->current_cmd_operation=="classify") RunClassify(&argv[next_arg],argc-next_arg,ENZYME_REG_FILE);
  else if (cmd_line->current_cmd_operation=="filter") RunFilter(&argv[next_arg],argc-next_arg,ENZYME_REG_FILE,STATS_FILE);
  else if (cmd_line->current_cmd_operation=="bin") RunBin(&argv[next_arg],argc-next_arg,BIN_SIZE,GENOME_REG_FILE,PRINT_MATRIX_FORMAT,SPLIT_MATRIX);
  else if (cmd_line->current_cmd_operation=="matrix") RunMatrix(&argv[next_arg],argc-next_arg,BIN_SIZE,MAX_DIST,ROTATE45,REF_REG_FILE,OUTPUT_PREFIX);
  else if (cmd_line->current_cmd_operation=="convert") RunConvert(COL_LABELS,SEPARATOR,NORM_CONSTANT,MIN_SCORE,MAX_DIST);
  else { cerr << "Unknown method '" << cmd_line->current_cmd_operation << "'!\n"; delete cmd_line; exit(1); }

  // clean up
  delete cmd_line;

  return 0;
}



