//
// Copyright (c) 2014 Aristotelis Tsirigos
// All rights reserved. This program and the accompanying materials are made available under the terms of the GNU General Public License v2.0 
// which accompanies this distribution, and is available at http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
//

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
#include "core.h"
#include "gtools-intervals.h"



using namespace std;


//---------------------------------------------------------------------------------//
// Types & constants                                                               //
//---------------------------------------------------------------------------------//

const string PROGRAM = "gtools-overlaps";
const long int BUFFER_SIZE = 10000;



  
//---------------------------------------------------------------------------------//
// Command-line options                                                            //
//---------------------------------------------------------------------------------//

bool HELP;
bool VERBOSE;
char *BIN_BITS;
bool IS_SORTED;
bool SORTED_BY_STRAND;
bool IGNORE_STRAND;
bool MERGE_LABELS;
char *LABEL_SEPARATOR;
bool PRINT_LABELS;
bool PRINT_REGIONS;
bool MATCH_GAPS;
long int MAX_LABEL_VALUE;
unsigned long int MIN_COUNT;
double MIN_RPKM;
double MIN_DENSITY;
char *OFFSET_OP;
bool OFFSET_FRACTION;
bool CENTER;
bool SUBSET_NONOVERLAPS;
bool OFFSET_SKIP_REF_GAPS;
char *QUERY_OP;
long int UPSTREAM_MAX_DISTANCE;
long int UPSTREAM_MIN_DISTANCE;
long int UPSTREAM_DIST;
long int DOWNSTREAM_DIST;
bool DISTANCE_FLAG;
long int PROXIMAL_DIST;
bool PRINT_HEADER;



//-------InitCmdLine-----------
//
CmdLineWithOperations *InitCmdLine(int argc, char *argv[], int *next_arg)
{
  // initialize
  CmdLineWithOperations *cmd_line = new CmdLineWithOperations(); 
  cmd_line->SetProgramName(PROGRAM,VERSION);

  // number of input files per operation
  map<string,int> min_files;

  // set operations
  cmd_line->AddOperation("annotate", "[OPTIONS] REFERENCE-REGION-FILE <TEST-REGION-FILE>", \
  "Annotates test regions according to reference regions.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operands: region, region-set\n\
  * Test region requirements: single-interval\n\
  * Reference region requirements: single-interval\n\
  * Region-set requirements: sorted if -S option is used"\
  );

  cmd_line->AddOperation("annotate2", "[OPTIONS] REFERENCE-REGION-FILE <TEST-REGION-FILE>", \
  "Annotates test regions according to reference regions (version 2).", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operands: region, region-set\n\
  * Test region requirements: single-interval\n\
  * Reference region requirements: single-interval\n\
  * Region-set requirements: sorted if -S option is used"\
  );

  cmd_line->AddOperation("bin", "[OPTIONS] REFERENCE-REGION-FILE <TEST-REGION-FILE>", \
  "Finds overlaps of interval pairs with reference regions.", \
  "* Test region file format: REG, GFF, BED, SAM (but only REG allows interchromosomal associations)\n\
  * Reference region file format: REG, GFF, BED, SAM\n\
  * Operands: interval-pairs, region-set\n\
  * Test region requirements: interval-pairs\n\
  * Reference region requirements: chromosome/strand-compatible, sorted and non-overlapping\n\
  * Test region-set requirements: none\n\
  * Reference region-set requirements: none"\
  );

  cmd_line->AddOperation("reduce", "[OPTIONS] REGION-FILE", \
  "Reduces number of input regions by removing regions covered by some higher-rank region.", \
  "* Region file format: REG, GFF, BED, SAM (but only REG allows interchromosomal associations)\n\
  * Reference region file format: REG, GFF, BED, SAM\n\
  * Operands: regions, region-set\n\
  * Region requirements: none\n\
  * Region-set requirements: none"\
  );

  cmd_line->AddOperation("count", "[OPTIONS] REFERENCE-REGION-FILE <TEST-REGION-FILE>", \
  "Counts the number of overlapping test regions per reference region.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operands: region, region-set\n\
  * Region requirements: chromosome/strand-compatible, sorted, non-overlapping\n\
  * Region-set requirements: sorted if -S option is used"\
  );

  cmd_line->AddOperation("coverage", "[OPTIONS] REFERENCE-REGION-FILE <TEST-REGION-FILE>", \
  "Calculates the depth coverage (i.e. the total number of overlapping nucleotides) per reference region.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operands: region, region-set\n\
  * Region requirements: chromosome/strand-compatible, sorted, non-overlapping\n\
  * Region-set requirements: sorted if -S option is used"\
  );

  cmd_line->AddOperation("density", "[OPTIONS] REFERENCE-REGION-FILE <TEST-REGION-FILE>", \
  "Computes the density (i.e. the coverage divided by the size of the reference region) of overlaps per reference region.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operands: region, region-set\n\
  * Region requirements: chromosome/strand-compatible, sorted, non-overlapping\n\
  * Region-set requirements: sorted if -S option is used"\
  );

  cmd_line->AddOperation("dist", "[OPTIONS] REFERENCE-REGION-FILE <TEST-REGION-FILE>", \
  "Computes the distance between a pair of intervals given breakpoints in reference file (e.g. restriction enzyme sites) [UNDER DEVELOPMENT].", \
  "* Test region file format: REG (because interchromosomal associations must be allowed)\n\
  * Reference region file format: REG, GFF, BED, SAM\n\
  * Operands: interval pair, region-set\n\
  * Test region requirements: none\n\
  * Reference region requirements: chromosome/strand-compatible, sorted and non-overlapping\n\
  * Test region-set requirements: none\n\
  * Reference region-set requirements: non-overlapping"\
  );

  cmd_line->AddOperation("intersect", "[OPTIONS] REFERENCE-REGION-FILE <TEST-REGION-FILE>", \
  "Computes the intersection between all pairs of test and reference regions. Results are grouped by test region.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operands: region, region-set\n\
  * Test region requirements: chromosome/strand-compatible, sorted, non-overlapping\n\
  * Reference region requirements: single-interval regions\n\
  * Region-set requirements: sorted if -S option is used"\
  );

  cmd_line->AddOperation("offset", "[OPTIONS] REFERENCE-REGION-FILE <TEST-REGION-FILE>", \
  "Computes the distances of test regions from their overlapping reference regions.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operands: region, region-set\n\
  * Test region requirements: single-interval\n\
  * Reference region requirements: chromosome/strand-compatible, sorted, non-overlapping\n\
  * Region-set requirements: sorted if -S option is used"\
  );

  cmd_line->AddOperation("overlap", "[OPTIONS] REFERENCE-REGION-FILE <TEST-REGION-FILE>", \
  "Finds the overlaps between all pairs of test and reference regions. Results are grouped by test region.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operands: region, region-set\n\
  * Region requirements: chromosome/strand-compatible, sorted, non-overlapping\n\
  * Region-set requirements: sorted if -S option is used"\
  );

  cmd_line->AddOperation("rpkm", "[OPTIONS] REFERENCE-REGION-FILE <TEST-REGION-FILE>", \
  "Computing reference region RPKM values.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operands: region, region-set\n\
  * Region requirements: chromosome/strand-compatible, sorted, non-overlapping\n\
  * Region-set requirements: sorted if -S option is used"\
  );

  cmd_line->AddOperation("subset", "[OPTIONS] REFERENCE-REGION-FILE <TEST-REGION-FILE>", \
  "Picks a subset of test regions depending on their overlap with reference regions. Results are grouped by test region.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operands: region, region-set\n\
  * Region requirements: chromosome/strand-compatible, sorted, non-overlapping\n\
  * Region-set requirements: sorted if -S option is used"\
  );

  if (argc<2) { cmd_line->OperationSummary("OPERATION [OPTIONS] REGION-FILE(S)","Performs overlap operations between a test and a reference set of genomic regions."); exit(1); }

  // current operation
  string op = argv[1];
  if (op[0]=='-') op = op.substr(1);		// ensure compatibility with previous version
  cmd_line->SetCurrentOperation(op);
  
  // common options
  cmd_line->AddOption("--help", &HELP, false, "help");
  cmd_line->AddOption("-h", &HELP, false, "help");
  cmd_line->AddOption("-v", &VERBOSE, false, "verbose mode");
  cmd_line->AddOption("-B", &BIN_BITS, "17,20,23,26", "number of shift-bits for each bin level");
  cmd_line->AddOption("-S", &IS_SORTED, false, "test and reference regions are sorted by chromosome and start position");
  cmd_line->AddOption("-s", &SORTED_BY_STRAND, false, "test and reference regions are also sorted by strand (-S must be set)");
  cmd_line->AddOption("-i", &IGNORE_STRAND, false, "ignore strand while finding overlaps");
   
  // Main options
  if (op=="annotate") {
    min_files["annotate"] = 1;
    cmd_line->AddOption("--query-op", &QUERY_OP, "overlap", "query operation for comparison with reference: {center|overlap}");
    cmd_line->AddOption("--upstream-max", &UPSTREAM_MAX_DISTANCE, 100000, "maximum allowed upstream region size");
    cmd_line->AddOption("--upstream-min", &UPSTREAM_MIN_DISTANCE, 10000, "minimum allowed upstream region size (subject to genomic bounds)");
    cmd_line->AddOption("--distance-flag", &DISTANCE_FLAG, false, "add proximal-distal indication");
    cmd_line->AddOption("--proximal-dist", &PROXIMAL_DIST, 1000, "define proximal distance (in nucleotides)");
    cmd_line->AddOption("--print-header", &PRINT_HEADER, false, "print header");
  }
  else if (op=="annotate2") {
    min_files["annotate2"] = 1;
    cmd_line->AddOption("--query-op", &QUERY_OP, "overlap", "query operation for comparison with reference: {center|overlap}");
    cmd_line->AddOption("--upstream-dist", &UPSTREAM_DIST, 100000, "maximum allowed upstream region size");
    cmd_line->AddOption("--downstream-dist", &DOWNSTREAM_DIST, 100000, "maximum allowed downstream region size");
    cmd_line->AddOption("--proximal-dist", &PROXIMAL_DIST, 1000, "define proximal distance (in nucleotides)");
    cmd_line->AddOption("--print-header", &PRINT_HEADER, false, "print header");
  }
  else if (op=="bin") {
    min_files["bin"] = 1;
    cmd_line->AddOption("--print-labels", &PRINT_LABELS, false, "print test region labels");
    cmd_line->AddOption("--print-regions", &PRINT_REGIONS, false, "print test regions");
  }
  else if (op=="reduce") {
    min_files["reduce"] = 1;
  }
  else if (op=="count") {
    min_files["count"] = 1;
    cmd_line->AddOption("-gaps", &MATCH_GAPS, false, "matching gaps between intervals are considered overlaps");
    cmd_line->AddOption("--max-label-value", &MAX_LABEL_VALUE, 1, "maximum region label value to be used");
    cmd_line->AddOption("-min", &MIN_COUNT, 0, "minimum count");
  }
  else if (op=="coverage") {
    min_files["coverage"] = 1;
    cmd_line->AddOption("-gaps", &MATCH_GAPS, false, "matching gaps between intervals are considered overlaps");
    cmd_line->AddOption("--max-label-value", &MAX_LABEL_VALUE, 1, "maximum region label value to be used");
    cmd_line->AddOption("-min", &MIN_COUNT, 0, "minimum coverage");
  }
  else if (op=="density") {
    min_files["density"] = 1;
    cmd_line->AddOption("-gaps", &MATCH_GAPS, false, "matching gaps between intervals are considered overlaps");
    cmd_line->AddOption("--max-label-value", &MAX_LABEL_VALUE, 1, "maximum region label value to be used");
    cmd_line->AddOption("-min", &MIN_DENSITY, 0.0, "minimum density");
  }
  else if (op=="dist") {
    min_files["dist"] = 1;
    cmd_line->AddOption("--print-labels", &PRINT_LABELS, false, "print test region labels");
    cmd_line->AddOption("--print-regions", &PRINT_REGIONS, false, "print test regions");
  }
  else if (op=="intersect") {
    min_files["intersect"] = 1;
    cmd_line->AddOption("-label", &MERGE_LABELS, false, "print query label for each match");
  }
  else if (op=="offset") {
    min_files["offset"] = 1;
    cmd_line->AddOption("-gaps", &MATCH_GAPS, false, "matching gaps between intervals are considered overlaps");
    cmd_line->AddOption("--skip-ref-gaps", &OFFSET_SKIP_REF_GAPS, false, "ignore gaps in reference regions when computing offsets");
    cmd_line->AddOption("-label", &PRINT_LABELS, false, "print test region labels");
    cmd_line->AddOption("-op", &OFFSET_OP, "5p", "reference point (1=start, 2=stop, 5p=5'-end, 3p=3'-end)");
    cmd_line->AddOption("-a", &OFFSET_FRACTION, false, "print distances as a fraction of total size");
    cmd_line->AddOption("-c", &CENTER, false, "print center of interval only");
  }
  else if (op=="overlap") {
    min_files["overlap"] = 1;
    cmd_line->AddOption("-gaps", &MATCH_GAPS, false, "matching gaps between intervals are considered overlaps");
    cmd_line->AddOption("-label", &MERGE_LABELS, false, "print query label for each match");
    cmd_line->AddOption("-t", &LABEL_SEPARATOR, ":", "label separator");
  }
  else if (op=="rpkm") {
    min_files["rpkm"] = 1;
    cmd_line->AddOption("-gaps", &MATCH_GAPS, false, "matching gaps between intervals are considered overlaps");
    cmd_line->AddOption("--max-label-value", &MAX_LABEL_VALUE, 1, "maximum region label value to be used");
    cmd_line->AddOption("-min", &MIN_RPKM, 0.0, "minimum RPKM");
  }
  else if (op=="subset") {
    min_files["subset"] = 1;
    cmd_line->AddOption("-gaps", &MATCH_GAPS, false, "matching gaps between intervals are considered overlaps");
    cmd_line->AddOption("-inv", &SUBSET_NONOVERLAPS, false, "print test regions that do *not* overlap with reference regions");
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




//-------PrintAnnotations-----------
//
void PrintAnnotations(GenomicRegion *qreg, GenomicRegion *ireg, bool ignore_strand, bool skip_ref_gaps, char *query_op, const char *offset_op, bool flag, long int proximal_dist)
{
  long int start_offset, stop_offset;
  qreg->I[0]->GetOffsetFrom(ireg,offset_op,ignore_strand,&start_offset,&stop_offset);
  double offset;
  if (strcmp(QUERY_OP,"center")==0) { offset = (double)(start_offset+stop_offset)/2; if (offset<0) return; }
  else if (strcmp(QUERY_OP,"overlap")==0) offset = (double)start_offset;
  else { fprintf(stderr, "Error [PrintAnnotations]: invalid value for query operation!\n"); exit(1); }
  //cout << "offsets: " << start_offset << ' ' << stop_offset << '\n';
  printf("%s\t", qreg->LABEL);
  qreg->I[0]->PrintInterval();
  printf("\t");
  printf("%lu\t", (unsigned long int)qreg->I[0]->GetSize());
  if ((strcmp(offset_op,"3p")==0)&&(flag==true)) printf("%s:", offset>=proximal_dist?"distal":"proximal");	
  else if ((strcmp(offset_op,"5p")==0)&&(flag==true)) printf("%s:downstream:", offset>proximal_dist?"distal":"proximal");
  printf("%s\t", ireg->LABEL);
  ireg->I[0]->PrintInterval();
  printf("\t");
  printf("%lu\t", (unsigned long int)ireg->I[0]->GetSize());
  printf("%ld\t", (long int)offset);
  printf("%f", offset/ireg->GetSize(skip_ref_gaps));
	printf("\n");
}






//---------CreateTerritory--------
//
GenomicRegionSet *CreateTerritory(GenomicRegionSet *RefRegSet, long int upstream_dist, long int downstream_dist)
{
	// constants
	const unsigned long int buffer_size = 10000;
	bool verbose = false;

	// create gene territory regions
	FILE *F = tmpfile();
	if (F==NULL) { cerr << "Error: cannot open temporary file!\n"; exit(1); }
	Progress PRG("Creating territory regions...",1);
  for (GenomicRegion *ireg=RefRegSet->Get(); ireg!=NULL; ireg=RefRegSet->Next()) {
    if (ireg->I.size()!=1) ireg->PrintError("single-interval reference regions are required for this operation!");
    char strand = ireg->I[0]->STRAND;
    long int start = ireg->I[0]->START;
    long int stop = ireg->I[0]->STOP;
    long int new_start = start - (strand=='+'?upstream_dist:downstream_dist);
    long int new_stop = stop + (strand=='+'?downstream_dist:upstream_dist);
  	fprintf(F, "%s\t%s %c %ld %ld\n", ireg->LABEL, ireg->I[0]->CHROMOSOME, strand, new_start, new_stop);
  	PRG.Check();
  }
  PRG.Done();
  RefRegSet->Reset();
  rewind(F);

  // return region set
  return new GenomicRegionSet(F,buffer_size,verbose,true,true);
}


















//---------------------------------------------------------------------------------//
// MAIN                                                                            //
//---------------------------------------------------------------------------------//
int main(int argc, char* argv[]) 
{
  // process command-line arguments
  int next_arg;
  CmdLineWithOperations *cmd_line = InitCmdLine(argc,argv,&next_arg); 
  _MESSAGES_ = VERBOSE;

  if (IS_SORTED&&SORTED_BY_STRAND&&IGNORE_STRAND) { fprintf(stderr, "[Error]: the input is sorted by chromosome/strand/start (i.e. -S and -s are set), therefore the overlap algorithm can only report strand-specific results (i.e. -i cannot be set)!\n"); exit(1); }

  //--------------------
  // annotate/unsorted
  //--------------------
  if (cmd_line->current_cmd_operation=="annotate") {
    // open region sets and index reference regions
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,true,true);
    GenomicRegionSetIndex RefIndex(&RefRegSet,BIN_BITS);
    bool skip_ref_gaps = false;
    bool match_gaps = true;		// NOTE: since we are only allowing single-interval regions, this is irrelevant

    // create upstream genomic annotator
    StringLIntMap *bounds = NULL;
    GenomicRegionSet *UpstreamRefRegSet = NULL;
    GenomicRegionSetIndex *UpstreamRefIndex = NULL;
    if (UPSTREAM_MAX_DISTANCE>0) {
    	UpstreamRefRegSet = CreateGenomicRegionSetAnnotator(&RefRegSet,bounds,IGNORE_STRAND,UPSTREAM_MAX_DISTANCE,UPSTREAM_MIN_DISTANCE,BIN_BITS);
    	UpstreamRefIndex = new GenomicRegionSetIndex(UpstreamRefRegSet,BIN_BITS);
    }

    // annotate
    Progress PRG("Annotating test regions...",1);
  	if (PRINT_HEADER) printf("TEST-LABEL\tTEST-LOCUS\tTEST-LOCUS-SIZE\tREF-LABEL\tREF-LOCUS\tREF-LOCUS-SIZE\tOFFSET\tNORMALIZED-OFFSET\n");
    for (GenomicRegion *qreg=TestRegSet.Get(); qreg!=NULL; qreg=TestRegSet.Next()) {
	  if (qreg->I.size()!=1) qreg->PrintError("single-interval test regions are required for this operation!");
	  GenomicInterval *qint = qreg->I[0];
	  for (GenomicRegion *ireg=RefIndex.GetOverlap(qint,match_gaps,IGNORE_STRAND); ireg!=NULL; ireg=RefIndex.NextOverlap(match_gaps,IGNORE_STRAND)) {
        if (ireg->I.size()!=1) ireg->PrintError("single-interval reference regions are required for this operation!");
        PrintAnnotations(qreg,ireg,IGNORE_STRAND,skip_ref_gaps,QUERY_OP,"5p",DISTANCE_FLAG,PROXIMAL_DIST);
	  }
	  if (UpstreamRefIndex) {
	    for (GenomicRegion *ireg=UpstreamRefIndex->GetOverlap(qint,match_gaps,IGNORE_STRAND); ireg!=NULL; ireg=UpstreamRefIndex->NextOverlap(match_gaps,IGNORE_STRAND)) {
          if (ireg->I.size()!=1) ireg->PrintError("single-interval reference regions are required for this operation!");
          PrintAnnotations(qreg,ireg,IGNORE_STRAND,skip_ref_gaps,QUERY_OP,"3p",DISTANCE_FLAG,PROXIMAL_DIST);
	    }
	  }
      PRG.Check();
    }
    PRG.Done();

    // cleanup
    if (UpstreamRefRegSet) delete UpstreamRefRegSet;
    if (UpstreamRefIndex) delete UpstreamRefIndex;

  }



  //---------------------------------------------
  // annotate/sorted
  //---------------------------------------------
  else if ((cmd_line->current_cmd_operation=="annotate")&&(IS_SORTED==true)) {
	  fprintf(stderr, "Error: operation annotate is not implemented for the -S option. Simply drop the -S and re-run!\n");
	  exit(1);
  }



  //--------------------
  // annotate2/unsorted
  //--------------------
  else if (cmd_line->current_cmd_operation=="annotate2") {
    if (IS_SORTED==true) fprintf(stderr, "Warning: ignoring the -S option.\n");
    
    // open region sets and index reference regions
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,true,true);
    bool match_gaps = true;		// NOTE: since we are only allowing single-interval regions, this is irrelevant

    // create upstream genomic annotator
    GenomicRegionSet *TerritoryRegSet = CreateTerritory(&RefRegSet,UPSTREAM_DIST,DOWNSTREAM_DIST);
    GenomicRegionSetIndex *TerritoryRegIndex = new GenomicRegionSetIndex(TerritoryRegSet,BIN_BITS);

    // annotate
    Progress PRG("Annotating test regions...",1);
  	if (PRINT_HEADER) printf("TEST-LABEL\tTEST-LOCUS\tTEST-LOCUS-SIZE\tREF-LABEL\tREGION\tDISTANCE-START\tDISTANCE-STOP\tREF-LOCUS\tREF-LOCUS-SIZE\n");
    for (GenomicRegion *qreg=TestRegSet.Get(); qreg!=NULL; qreg=TestRegSet.Next()) {
      if (qreg->I.size()!=1) qreg->PrintError("single-interval test regions are required for this operation!");
      GenomicInterval *qint = qreg->I[0];
      GenomicRegion *rreg = TerritoryRegIndex->GetOverlap(qint,match_gaps,IGNORE_STRAND);
      if (rreg==NULL) {
        printf("%s", qreg->LABEL);
        printf("\t");
        printf("%s:%ld-%ld", qint->CHROMOSOME, qint->START, qint->STOP);
        printf("\t");
        printf("%lu", (unsigned long int)qint->GetSize());
        printf("\t");
        printf("NA");
        printf("\t");
        printf("intergenic");
        printf("\t");
        printf("NA");
        printf("\t");
        printf("NA");
        printf("\t");
        printf("NA");
        printf("\t");
        printf("NA");
	      printf("\n");
      }
	    for ( ; rreg!=NULL; rreg=TerritoryRegIndex->NextOverlap(match_gaps,IGNORE_STRAND)) {
        if (rreg->I.size()!=1) rreg->PrintError("single-interval reference regions are required for this operation!");
        GenomicInterval *rint = rreg->I[0];
        long tss = rint->STRAND=='+' ? (rint->START+UPSTREAM_DIST) : (rint->STOP-UPSTREAM_DIST);
        long tes = rint->STRAND=='+' ? (rint->STOP-DOWNSTREAM_DIST) : (rint->START+DOWNSTREAM_DIST);
        GenomicInterval zint(rint->CHROMOSOME,rint->STRAND,rint->STRAND=='+'?tss:tes,rint->STRAND=='+'?tes:tss);                // re-create original reference region
        long dist_tss_start = rint->STRAND=='+' ? (qint->START-tss) : (tss - qint->STOP);
        long dist_tss_stop = rint->STRAND=='+' ? (qint->STOP-tss) : (tss - qint->START);
        printf("%s", qreg->LABEL);
        printf("\t");
        printf("%s:%ld-%ld", qint->CHROMOSOME, qint->START, qint->STOP);
        printf("\t");
        printf("%lu", (unsigned long int)qint->GetSize());
        printf("\t");
        printf("%s", rreg->LABEL);
        printf("\t");
        bool proximal = min(abs(dist_tss_start),abs(dist_tss_stop))<=PROXIMAL_DIST;
        bool upstream = min(dist_tss_start,dist_tss_stop)<-PROXIMAL_DIST;
        bool downstream = max(dist_tss_start,dist_tss_stop)>PROXIMAL_DIST;
        if (proximal) printf("proximal");
        else if (upstream&&downstream) printf("proximal");
        else if (upstream) printf("upstream");
        else if (downstream) printf("downstream"); 
        else printf("unclassified");
        printf("\t");
        printf("%ld", dist_tss_start);
        printf("\t");
        printf("%ld", dist_tss_stop);
        printf("\t");
        zint.PrintInterval();
        printf("\t");
        printf("%lu", (unsigned long int)zint.GetSize());
	      printf("\n");
	    }
      PRG.Check();
    }
    PRG.Done();

    // cleanup
    if (TerritoryRegSet) delete TerritoryRegSet;
    if (TerritoryRegIndex) delete TerritoryRegIndex;

  }



  //--------------------
  // bin/unsorted
  //--------------------
  else if (cmd_line->current_cmd_operation=="bin") {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,true,true);

    // classify into bins
    GenomicRegionSetIndex RefIndex(&RefRegSet,BIN_BITS);
    Progress PRG("Classifying interval pairs into reference bins...",1);
	  bool match_gaps = false;
    for (GenomicRegion *qreg=TestRegSet.Get(); qreg!=NULL; qreg=TestRegSet.Next()) {
	  if (qreg->I.size()!=2) qreg->PrintError("interval pairs are required for this operation!");
	  GenomicInterval *qint1 = qreg->I[0];
	  GenomicInterval *qint2 = qreg->I[1];
	  typedef list<GenomicRegion*> GenomicRegionList;
	  GenomicRegionList regList1, regList2;
    for (GenomicRegion *ireg1=RefIndex.GetOverlap(qint1,match_gaps,IGNORE_STRAND); ireg1!=NULL; ireg1=RefIndex.NextOverlap(match_gaps,IGNORE_STRAND)) regList1.push_back(ireg1);
    for (GenomicRegion *ireg2=RefIndex.GetOverlap(qint2,match_gaps,IGNORE_STRAND); ireg2!=NULL; ireg2=RefIndex.NextOverlap(match_gaps,IGNORE_STRAND)) regList2.push_back(ireg2);
	  for (GenomicRegionList::iterator x=regList1.begin(); x!=regList1.end(); x++)
	    for (GenomicRegionList::iterator y=regList2.begin(); y!=regList2.end(); y++) {
        printf("%s\t", qreg->LABEL); 
	      (*x)->I[0]->PrintInterval(); // NOTE: print entire region
		    printf(" "); 
		    (*y)->I[0]->PrintInterval(); // NOTE: print entire region
		    if (PRINT_LABELS) printf("\t%s\t%s", (*x)->LABEL, (*y)->LABEL); 
		    if (PRINT_REGIONS) { printf("\t"); qint1->PrintInterval(); printf(" "); qint2->PrintInterval(); }
		    printf("\n");
	    }
      PRG.Check();
    }
    PRG.Done();
  }


  //--------------------
  // reduce/unsorted
  //--------------------
  else if (cmd_line->current_cmd_operation=="reduce") {
    // open region sets
    char *REG_FILE = argv[next_arg];
    GenomicRegionSet RegSet(REG_FILE,BUFFER_SIZE,VERBOSE,true,true);

    // create a new region set that contains all intervals 
    Progress PRG1("Preprocessing input regions...",RegSet.n_regions);
    long int n_intervals = 0;
    for (long int k=0; k<RegSet.n_regions; k++) n_intervals += RegSet.R[k]->I.size();
    GenomicRegionSet IntervalSet(n_intervals,RegSet.format,VERBOSE);
    for (long int k=0,z=0; k<RegSet.n_regions; k++) {
      GenomicRegion *r = RegSet.R[k];
      for (size_t i=0; i<r->I.size(); i++) {
        IntervalSet.R[z] = new GenomicRegion(r->LABEL,r->I[i]);
        IntervalSet.R[z]->n_line = k;
        z++;
      }
      PRG1.Check();
    }
    PRG1.Done();
    
    // reduce input regions by removing regions that are "covered" by some higher-rank region
    GenomicRegionSetIndex IntervalIndex(&IntervalSet,BIN_BITS);
    Progress PRG2("Reducing regions...",1);
	  bool match_gaps = false;
    for (long int kquery=0; kquery<RegSet.n_regions; kquery++) {
      GenomicRegion *qreg = RegSet.R[kquery];
      bool redundant = false;
      for (size_t q=0; (redundant==false)&&(q<qreg->I.size()); q++) {
        GenomicInterval *qint = qreg->I[q];
        for (GenomicRegion *ireg=IntervalIndex.GetOverlap(qint,match_gaps,IGNORE_STRAND); (redundant==false)&&(ireg!=NULL); ireg=IntervalIndex.NextOverlap(match_gaps,IGNORE_STRAND)) {
          long int kmatch = ireg->n_line;
          if ((kquery>kmatch)&&(RegSet.R[kquery]->IsCoveredBy(RegSet.R[kmatch],IGNORE_STRAND)==true)) redundant = true;
        }
      } 
      if (redundant==false) RegSet.R[kquery]->Print(); 
      PRG2.Check();
    }
    PRG2.Done();
  }


  //--------------------
  // count
  //--------------------
  else if (cmd_line->current_cmd_operation=="count") {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,true,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);

    GenomicRegionSetOverlaps *overlaps;
    if (IS_SORTED) overlaps = new SortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,SORTED_BY_STRAND);
    else overlaps = new UnsortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,BIN_BITS);
    unsigned long int *hits = overlaps->CountIndexOverlaps(MATCH_GAPS,IGNORE_STRAND,MAX_LABEL_VALUE); 
    Progress PRG("Printing counts...",RefRegSet.n_regions);
    for (long int k=0; k<RefRegSet.n_regions; k++) {
      if (hits[k]>=MIN_COUNT) {
        printf("%s", RefRegSet.R[k]->LABEL); printf("\t");
        printf("%lu", hits[k]);
        printf("\n");
      }
      PRG.Check();
    }
    PRG.Done();
    delete hits;
    delete overlaps;
  }



  //--------------------
  // coverage
  //--------------------
  else if (cmd_line->current_cmd_operation=="coverage") {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,true,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);

    GenomicRegionSetOverlaps *overlaps;
    if (IS_SORTED) overlaps = new SortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,SORTED_BY_STRAND);
    else overlaps = new UnsortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,BIN_BITS);
    unsigned long int *coverage = overlaps->CalcIndexCoverage(MATCH_GAPS,IGNORE_STRAND,MAX_LABEL_VALUE); 
    Progress PRG("Printing coverages...",RefRegSet.n_regions);
    for (long int k=0; k<RefRegSet.n_regions; k++) {
      if (coverage[k]>=MIN_COUNT) {
        printf("%s", RefRegSet.R[k]->LABEL); printf("\t");
        printf("%lu", coverage[k]); printf("\n");
      }
      PRG.Check();
    }
    PRG.Done();
    delete coverage;
    delete overlaps;
  }



  //--------------------
  // density
  //--------------------
  else if (cmd_line->current_cmd_operation=="density") {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,true,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);

    // process overlaps
    GenomicRegionSetOverlaps *overlaps;
    if (IS_SORTED) overlaps = new SortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,SORTED_BY_STRAND);
    else overlaps = new UnsortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,BIN_BITS);
    unsigned long int *coverage = overlaps->CalcIndexCoverage(MATCH_GAPS,IGNORE_STRAND,MAX_LABEL_VALUE); 
    Progress PRG("Printing densities...",RefRegSet.n_regions);
    for (long int k=0; k<RefRegSet.n_regions; k++) {
      GenomicRegion *qreg = RefRegSet.R[k];
      long int qreg_size = qreg->GetSize(!MATCH_GAPS);
      double density = (double)coverage[k]/qreg_size;
      if (density>=MIN_DENSITY) printf("%s\t%.4e\n", qreg->LABEL, density);
      PRG.Check();
    }
    PRG.Done();
    delete coverage;
    delete overlaps;
  }


  //--------------------
  // dist/unsorted
  //--------------------
  else if (cmd_line->current_cmd_operation=="dist") {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);
    if (TestRegSet.format!="REG") { fprintf(stderr, "Error: test region set should be in REG format for this operation!\n"); exit(1); }
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,true,true);

    // compute distances
    GenomicRegionSetIndex RefIndex(&RefRegSet,BIN_BITS);
    Progress PRG("Processing interval pairs...",1);
    for (GenomicRegion *qreg=TestRegSet.Get(); qreg!=NULL; qreg=TestRegSet.Next()) {
	  if (qreg->I.size()!=2) qreg->PrintError("interval pairs are required for this operation!");
	  const char *op[] = { "3p", "5p" };
	  bool matched = false;
	  long int offsets[4];
	  for (int t=0,k=0; t<=1; t++) {
	    for (int q=0; q<=1; q++,k++) {
	      GenomicInterval interval = qreg->I[q];
	      interval.ModifyPos(op[t],0);
          for (GenomicRegion *ireg=RefIndex.GetMatch(&interval); ireg!=NULL; ireg=RefIndex.NextMatch()) {
            long int start_offset, stop_offset;
            interval.GetOffsetFrom(ireg->I.front(),op[1-t],true,&start_offset,&stop_offset);
            offsets[k] = start_offset;
			matched = true;
		  }
        }
	  }
	  if (matched) {
        if (PRINT_LABELS) printf("%s\t", qreg->LABEL);
	    if (PRINT_REGIONS) { qreg->I[0]->PrintInterval(); printf(" "); qreg->I[1]->PrintInterval(); printf("\t"); }
	    for (int k=0; k<4; k++) printf("%ld%c", offsets[k], k==3?'\t':' ');
		//if (qreg->I[0]->STRAND==qreg->I[1]->STRAND) printf("%ld", min(offsets[0]+offsets[1],offsets[2]+offsets[3]));
		//else printf("%ld", min(offsets[0]+offsets[3],offsets[1]+offsets[2]));
		long int d1 = min(offsets[0]+offsets[1],offsets[2]+offsets[3]);
		long int d2 = min(offsets[0]+offsets[3],offsets[1]+offsets[2]);
		printf("%ld", min(d1,d2));
		printf("\t");
		if (strcmp(qreg->I[0]->CHROMOSOME,qreg->I[1]->CHROMOSOME)==0) printf("%ld", max(qreg->I[0]->STOP,qreg->I[1]->STOP)-min(qreg->I[0]->START,qreg->I[1]->START));
		else printf("Inf");
        printf("\n");
	  }
      PRG.Check();
    }
    PRG.Done();
  }


  //---------------------------------------------
  // offset/sorted/no-subtract-gaps
  //---------------------------------------------
  else if ((cmd_line->current_cmd_operation=="offset")&&(IS_SORTED==true)&&(OFFSET_SKIP_REF_GAPS==false)) {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);

    // compute offsets
    SortedGenomicRegionSetOverlaps Overlaps(&RefRegSet,&TestRegSet,SORTED_BY_STRAND);
    Progress PRG("Processing queries...",1);
    for (GenomicRegion *qreg=Overlaps.GetQuery(); Overlaps.Done()==false; qreg=Overlaps.NextQuery()) {
      size_t Qsize = qreg->GetSize(OFFSET_SKIP_REF_GAPS);
      for (GenomicRegion *ireg=Overlaps.GetOverlap(MATCH_GAPS,IGNORE_STRAND); ireg!=NULL; ireg=Overlaps.NextOverlap(MATCH_GAPS,IGNORE_STRAND)) {
        if (ireg->I.size()>1) ireg->PrintError("multi-interval test regions are not allowed for this operation!"); // NOTE: we may want to remove this restriction in the future
        cout << qreg->LABEL << '\t';
        if (PRINT_LABELS) cout << ireg->LABEL << ' ';
        long int start_offset, stop_offset;
        GenomicInterval iint(ireg->I.front()->CHROMOSOME,ireg->I.front()->STRAND,ireg->I.front()->START,ireg->I.back()->STOP);
        iint.GetOffsetFrom(qreg,OFFSET_OP,IGNORE_STRAND,&start_offset,&stop_offset);
        if (start_offset>stop_offset) { fprintf(stderr, "Error: start offset is greater than stop offset (this must be a bug)!\n"); exit(1); }
        if (CENTER) {
          if (OFFSET_FRACTION) printf("%f", ((float)start_offset/Qsize+(float)stop_offset/Qsize)/2);
          else printf("%ld", (start_offset+stop_offset)/2);
        }
        else {
          if (OFFSET_FRACTION) printf("%f %f", (float)start_offset/Qsize, (float)stop_offset/Qsize);
          else printf("%ld %ld", start_offset, stop_offset);
        }
        cout << '\n';
      }
      PRG.Check();
    }
    PRG.Done();
  }



  //---------------------------------------------
  // offset/sorted/subtract-gaps
  //---------------------------------------------
  else if ((cmd_line->current_cmd_operation=="offset")&&(IS_SORTED==true)&&(OFFSET_SKIP_REF_GAPS==true)) {
	  fprintf(stderr, "Error: option --subtract-gaps is not implemented for the -S option. Simply drop the -S and re-run!\n");
	  exit(1);
  }



  //---------------------------------------------
  // offset/unsorted/no-subtract-gaps
  //---------------------------------------------
  else if ((cmd_line->current_cmd_operation=="offset")&&(IS_SORTED==false)&&(OFFSET_SKIP_REF_GAPS==false)) {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,true,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);

    // compute offsets
    UnsortedGenomicRegionSetOverlaps Overlaps(&TestRegSet,&RefRegSet,BIN_BITS);
    Progress PRG("Processing queries...",1);
    for (GenomicRegion *qreg=Overlaps.GetQuery(); Overlaps.Done()==false; qreg=Overlaps.NextQuery()) {
      for (GenomicRegion *ireg=Overlaps.GetOverlap(MATCH_GAPS,IGNORE_STRAND); ireg!=NULL; ireg=Overlaps.NextOverlap(MATCH_GAPS,IGNORE_STRAND)) {
        if (qreg->I.size()>1) qreg->PrintError("multi-interval test regions are not allowed for this operation!"); // NOTE: we may want to remove this restriction in the future
        size_t Isize = ireg->GetSize(OFFSET_SKIP_REF_GAPS);
        printf("%s\t", ireg->LABEL);
        if (PRINT_LABELS) printf("%s ", qreg->LABEL);
        long int start_offset, stop_offset;
        GenomicInterval qint(qreg->I.front()->CHROMOSOME,qreg->I.front()->STRAND,qreg->I.front()->START,qreg->I.back()->STOP);
        qint.GetOffsetFrom(ireg,OFFSET_OP,IGNORE_STRAND,&start_offset,&stop_offset);
        if (start_offset>stop_offset) { fprintf(stderr, "Error: start offset is greater than stop offset (this must be a bug)!\n"); exit(1); }
        if (CENTER) {
          if (OFFSET_FRACTION) printf("%f", ((float)start_offset/Isize+(float)stop_offset/Isize)/2);
          else printf("%ld", (start_offset+stop_offset)/2);
        }
        else {
          if (OFFSET_FRACTION) printf("%f %f", (float)start_offset/Isize, (float)stop_offset/Isize);
          else printf("%ld %ld", start_offset, stop_offset);
        }
        printf("\n");
      }
      PRG.Check();
    }
    PRG.Done();
  }


  //---------------------------------------------
  // offset/unsorted/subtract-gaps
  //---------------------------------------------
  else if ((cmd_line->current_cmd_operation=="offset")&&(IS_SORTED==false)&&(OFFSET_SKIP_REF_GAPS==true)) {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,true,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);

    // compute offsets
    UnsortedGenomicRegionSetOverlaps Overlaps(&TestRegSet,&RefRegSet,BIN_BITS);
    Progress PRG("Processing queries...",1);
    for (GenomicRegion *qreg=Overlaps.GetQuery(); Overlaps.Done()==false; qreg=Overlaps.NextQuery()) {
      for (GenomicRegion *ireg=Overlaps.GetOverlap(MATCH_GAPS,IGNORE_STRAND); ireg!=NULL; ireg=Overlaps.NextOverlap(MATCH_GAPS,IGNORE_STRAND)) {
    	  OffsetList *offsets = CalcOffsetsWithoutGaps(qreg,ireg,OFFSET_OP,IGNORE_STRAND);
          if ((offsets!=NULL)&&(offsets->size()>0)) {
              printf("%s\t", ireg->LABEL);
              if (PRINT_LABELS) printf("%s ", qreg->LABEL);
              size_t Isize = ireg->GetSize(OFFSET_SKIP_REF_GAPS);
              for (OffsetList::iterator x=offsets->begin(); x!=offsets->end(); x++) {
                  if (CENTER) {
                      if (OFFSET_FRACTION) printf("%f", ((float)x->first/Isize+(float)x->second/Isize)/2);
                      else printf("%ld", (x->first+x->second)/2);
                  }
                  else {
                      if (OFFSET_FRACTION) printf("%f %f", (float)x->first/Isize, (float)x->second/Isize);
                      else printf("%ld %ld", x->first, x->second);
                  }
              }
              printf("\n");
          }
          if (offsets!=NULL) delete offsets;
      }
      PRG.Check();
    }
    PRG.Done();
  }


  //--------------------
  // intersect
  //--------------------
  else if (cmd_line->current_cmd_operation=="intersect") {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,IS_SORTED?false:true,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,false);

    // process
    GenomicRegionSetOverlaps *overlaps;
    if (IS_SORTED) overlaps = new SortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,SORTED_BY_STRAND);
    else overlaps = new UnsortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,BIN_BITS);
    Progress PRG("Processing queries...",1);
    for (GenomicRegion *qreg=overlaps->GetQuery(); overlaps->Done()==false; qreg=overlaps->NextQuery()) {   
      GenomicRegion *ireg = overlaps->GetOverlap(false,IGNORE_STRAND);
      if (ireg!=NULL) {
        do { 
          qreg->PrintConstrained(ireg,MERGE_LABELS);
          ireg = overlaps->NextOverlap(false,IGNORE_STRAND);
        } while (ireg!=NULL);
      }
      PRG.Check();
    }
    PRG.Done();
    delete overlaps;
  }


  //--------------------
  // overlap
  //--------------------
  else if (cmd_line->current_cmd_operation=="overlap") {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,IS_SORTED?false:true,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,false);

    // process
    GenomicRegionSetOverlaps *overlaps;
    if (IS_SORTED) overlaps = new SortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,SORTED_BY_STRAND);
    else overlaps = new UnsortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,BIN_BITS);
    Progress PRG("Computing query overlaps...",1);
    for (GenomicRegion *qreg=overlaps->GetQuery(); overlaps->Done()==false; qreg=overlaps->NextQuery()) {
      GenomicRegion *ireg = overlaps->GetOverlap(MATCH_GAPS,IGNORE_STRAND);
      if (ireg!=NULL) {
        do { 
          if (MERGE_LABELS) {
            char *old_label = qreg->LABEL;
            char *new_label = new char[strlen(qreg->LABEL)+1+strlen(ireg->LABEL)+1];
            sprintf(new_label, "%s%s%s", qreg->LABEL, LABEL_SEPARATOR, ireg->LABEL);
            qreg->LABEL = new_label;
            qreg->Print();
            qreg->LABEL = old_label;
            delete new_label;
          }
          else qreg->Print();
          ireg = overlaps->NextOverlap(MATCH_GAPS,IGNORE_STRAND);
        } while (ireg!=NULL);
      }
      PRG.Check();
    }
    PRG.Done();
    delete overlaps;
  }

  

  //--------------------
  // rpkm
  //--------------------
  else if (cmd_line->current_cmd_operation=="rpkm") {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,true,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,true);

    GenomicRegionSetOverlaps *overlaps;
    if (IS_SORTED) overlaps = new SortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,SORTED_BY_STRAND);
    else overlaps = new UnsortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,BIN_BITS);
    unsigned long int *hits = overlaps->CountIndexOverlaps(MATCH_GAPS,IGNORE_STRAND,MAX_LABEL_VALUE); 
    unsigned long int nreads = 0;
    for (long int k=0; k<RefRegSet.n_regions; k++) nreads += hits[k]; 
    if (VERBOSE) fprintf(stderr, "* %lu reads overlap reference regions.\n", nreads);
    double mreads = (double)nreads/1000000;
    Progress PRG("Printing counts...",RefRegSet.n_regions);
    for (long int k=0; k<RefRegSet.n_regions; k++) {
      if (hits[k]>=MIN_COUNT) {
        printf("%s", RefRegSet.R[k]->LABEL); printf("\t");
        long int eff_len = RefRegSet.R[k]->GetSize(!MATCH_GAPS);
        double rpkm = eff_len<=0?0.0/0.0:(double)1000*hits[k]/eff_len/mreads;
        printf("%.4e", rpkm);
        printf("\n");
      }
      PRG.Check();
    }
    PRG.Done();
    delete hits;
    delete overlaps;
  }



  //--------------------
  // subset
  //--------------------
  else if (cmd_line->current_cmd_operation=="subset") {
    // open region sets
    char *REF_REG_FILE = argv[next_arg];
    char *TEST_REG_FILE = next_arg+1==argc ? NULL : argv[next_arg+1];
    GenomicRegionSet RefRegSet(REF_REG_FILE,BUFFER_SIZE,VERBOSE,IS_SORTED?false:true,true);
    GenomicRegionSet TestRegSet(TEST_REG_FILE,BUFFER_SIZE,VERBOSE,false,false);

    // process
    GenomicRegionSetOverlaps *overlaps;
    if (IS_SORTED) overlaps = new SortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,SORTED_BY_STRAND);
    else overlaps = new UnsortedGenomicRegionSetOverlaps(&TestRegSet,&RefRegSet,BIN_BITS);
    Progress PRG("Creating query subset...",1);
    for (GenomicRegion *qreg=overlaps->GetQuery(); (SUBSET_NONOVERLAPS&&(qreg!=NULL))||(overlaps->Done()==false); qreg=overlaps->NextQuery()) {   
      if ((overlaps->GetOverlap(MATCH_GAPS,IGNORE_STRAND)==NULL)==SUBSET_NONOVERLAPS) qreg->Print(); 
      PRG.Check();
    }
    PRG.Done();
    delete overlaps;
  }



  // clean up
  delete cmd_line;
  
  return 0;
}



