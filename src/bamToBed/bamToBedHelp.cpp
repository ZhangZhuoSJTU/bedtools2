/*****************************************************************************
  bamToBedHelp.cpp

  Help text for bedtools bamtobed, split out of bamToBed.cpp
******************************************************************************/
#include "version.h"
#include <iostream>

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools bamtobed"

void bamtobed_help(void) {
    
    cerr << "\nTool:    bedtools bamtobed (aka bamToBed)" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Converts BAM alignments to BED6 or BEDPE format." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bam> " << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-bedpe\t"      << "Write BEDPE format." << endl;
    cerr                      << "\t\t- Requires BAM to be grouped or sorted by query." << endl << endl;

    cerr << "\t-mate1\t"      << "When writing BEDPE (-bedpe) format, " << endl;
    cerr                      << "\t\talways report mate one as the first BEDPE \"block\"." << endl << endl;
  
    cerr << "\t-bed12\t"      << "Write \"blocked\" BED format (aka \"BED12\"). Forces -split." << endl << endl;
    cerr                      << "\t\thttp://genome-test.cse.ucsc.edu/FAQ/FAQformat#format1" << endl << endl;
  
    cerr << "\t-split\t"      << "Report \"split\" BAM alignments as separate BED entries." << endl ;
    cerr                      << "\t\tSplits only on N CIGAR operations." << endl << endl;

    cerr << "\t-splitD\t"   << "Split alignments based on N and D CIGAR operators." << endl;
    cerr                      << "\t\tForces -split." << endl << endl;
  
    cerr << "\t-ed\t"         << "Use BAM edit distance (NM tag) for BED score." << endl;
    cerr                      << "\t\t- Default for BED is to use mapping quality." << endl;
    cerr                      << "\t\t- Default for BEDPE is to use the minimum of" << endl;
    cerr                      << "\t\t  the two mapping qualities for the pair." << endl;
    cerr                      << "\t\t- When -ed is used with -bedpe, the total edit" << endl;
    cerr                      << "\t\t  distance from the two mates is reported." << endl << endl;
  
    cerr << "\t-tag\t"        << "Use other NUMERIC BAM alignment tag for BED score." << endl;
    cerr                      << "\t\t- Default for BED is to use mapping quality." << endl;
    cerr                      << "\t\t  Disallowed with BEDPE output." << endl << endl;
  
    cerr << "\t-color\t"      << "An R,G,B string for the color used with BED12 format." << endl;
    cerr                      << "\t\tDefault is (255,0,0)." << endl << endl;
  
    cerr << "\t-cigar\t"      << "Add the CIGAR string to the BED entry as a 7th column." << endl << endl;


    // end the program here
    exit(1);
}
