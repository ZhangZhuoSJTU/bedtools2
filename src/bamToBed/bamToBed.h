#ifndef BAMTOBED_H
#define BAMTOBED_H

#include "api/BamReader.h"
#include "api/BamAux.h"
#include "BlockedIntervals.h"
#include "bedFile.h"
using namespace BamTools;

#include <vector>
#include <algorithm>    // for swap()
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>

using namespace std;

// function declarations
void bamtobed_help(void);

void ConvertBamToBed(const string &bamFile, bool useEditDistance, 
                     const string &bamTag, bool writeBed12, 
                     bool obeySplits, bool splitOnDeletions, 
                     const string &color, bool useCigar,
                     bool useNovoalign, bool useBWA);
                     
void ConvertBamToBedpe(const string &bamFile, 
                       const bool &useEditDistance, bool mate1First);

string PrintTag(const BamAlignment &bam, const string &tag);

void PrintBed(const BamAlignment &bam, const RefVector &refs, 
              bool useEditDistance, const string &bamTag, 
              bool obeySplits, bool splitOnDeletions, 
              bool useCigar, bool useNovoalign, 
              bool useBWA);
              
void PrintBed12(const BamAlignment &bam, const RefVector &refs, 
                bool useEditDistance, const string &bamTag,
                bool obeySplits, bool splitOnDeletions, 
                string color = "255,0,0");

void PrintBedPE(const BamAlignment &bam1, const BamAlignment &bam2,
                const RefVector &refs, bool useEditDistance, bool mate1First);

void ParseCigarBed12(const vector<CigarOp> &cigar, vector<int> &blockStarts,
                     vector<int> &blockEnds, int &alignmentEnd);

string BuildCigarString(const vector<CigarOp> &cigar);

bool bamtobed_IsCorrectMappingForBEDPE (const BamAlignment &bam);

#endif /* BAMTOBED_H */
