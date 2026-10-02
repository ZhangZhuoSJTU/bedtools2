/*****************************************************************************
  bamToBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "api/BamReader.h"
#include "api/BamAux.h"
#include "BlockedIntervals.h"
#include "bedFile.h"
#include "version.h"
using namespace BamTools;

#include <vector>
#include <algorithm>    // for swap()
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>

using namespace std;


// define our program name
#define PROGRAM_NAME "bedtools bamtobed"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)


#include "bamToBed.h"

int bamtobed_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bamFile = "stdin";
    string color   = "255,0,0";
    string tag     = "";

    bool haveBam           = true;
    bool haveColor         = false;
    bool haveOtherTag      = false;
    bool writeBedPE        = false;
    bool writeBed12        = false;
    bool useEditDistance   = false;
    bool useAlignmentScore = false;
    bool useCigar          = false;
    bool obeySplits        = false;
    bool splitOnDeletions  = false;
    bool mate1First        = false;
    bool useNovoalign      = false;  // custom for Quinlan/Hall research
    bool useBWA            = false;  // custom for Quinlan/Hall research

    // check to see if we should print out some help

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) bamtobed_help();

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                bamFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-bedpe", 6, parameterLength)) {
            writeBedPE = true;
        }
        else if(PARAMETER_CHECK("-bed12", 6, parameterLength)) {
            writeBed12 = true;
        }
        else if(PARAMETER_CHECK("-split", 6, parameterLength)) {
            obeySplits = true;
        }
        else if(PARAMETER_CHECK("-splitD", 7, parameterLength)) {
            obeySplits = true;
            splitOnDeletions = true;
        }
        else if(PARAMETER_CHECK("-ed", 3, parameterLength)) {
            useEditDistance = true;
            tag = "NM";
        }
        else if(PARAMETER_CHECK("-cigar", 6, parameterLength)) {
            useCigar = true;
        }
        else if(PARAMETER_CHECK("-as", 3, parameterLength)) {
            useAlignmentScore = true;
        }
        else if(PARAMETER_CHECK("-novo", 5, parameterLength)) {
            useNovoalign = true;
        }
        else if(PARAMETER_CHECK("-bwa", 4, parameterLength)) {
            useBWA = true;
        }
        else if(PARAMETER_CHECK("-mate1", 6, parameterLength)) {
            mate1First = true;
        }
        else if(PARAMETER_CHECK("-color", 6, parameterLength)) {
            if ((i+1) < argc) {
                haveColor = true;
                color = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-tag", 4, parameterLength)) {
            if ((i+1) < argc) {
                haveOtherTag = true;
                tag = argv[i + 1];
                i++;
            }
        }
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have an input files
    if (haveBam == false) {
        cerr << endl << "*****" << endl 
             << "*****ERROR: Need -i (BAM) file. " 
             << endl << "*****" << endl;
        showHelp = true;
    }
    if (haveColor == true && writeBed12 == false) {
        cerr << endl << "*****" << endl 
             << "*****ERROR: Cannot use color without BED12. " 
             << endl << "*****" << endl;
        showHelp = true;
    }
    if (useEditDistance == true && obeySplits == true) {
        cerr << endl << "*****" << endl 
            << "*****ERROR: Cannot use -ed with -splits. Edit distances cannot be computed for each \'chunk\'." 
            << endl << "*****" << endl;
        showHelp = true;
    }
    if (useEditDistance == true && useCigar == true) {
        cerr << endl << "*****" << endl 
            << "*****ERROR: Cannot use -cigar with -splits.  Not yet supported." << endl 
            << "*****" << endl;
        showHelp = true;
    }
    if (useEditDistance == true && haveOtherTag == true) {
        cerr << endl << "*****" << endl 
            << "*****ERROR: Cannot use -ed with -tag.  Choose one or the other." << endl 
            << "*****" << endl;
        showHelp = true;
    }
    if (writeBedPE == true && haveOtherTag == true) {
        cerr << endl << "*****" << endl 
            << "*****ERROR: Cannot use -tag with -bedpe." << endl 
            << "*****" << endl;
        showHelp = true;
    }
    if (writeBedPE == false && mate1First == true) {
        cerr << endl << "*****" << endl 
            << "*****ERROR: Must use -mate1 with -bedpe." << endl 
            << "*****" << endl;
        showHelp = true;
    }
    // if there are no problems, let's convert BAM to BED or BEDPE
    if (!showHelp) {
        if (writeBedPE == false)
            ConvertBamToBed(bamFile, useEditDistance, 
                            tag, writeBed12, 
                            obeySplits, splitOnDeletions, 
                            color, useCigar,
                            useNovoalign, useBWA);
        else
            ConvertBamToBedpe(bamFile, useEditDistance, mate1First);
    }
    else {
        bamtobed_help();
    }
    return 0;
}
