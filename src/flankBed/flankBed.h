/*****************************************************************************
  flankBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/

#include "bedFile.h"
#include "GenomeFile.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <cstdlib>
#include <ctime>
using namespace std;


//************************************************
// Named constants (replace hard-coded literals).
//************************************************
namespace FlankBedConstants {
    // Sentinel returned by GenomeFile::getChromSize when a chromosome is
    // absent from the genome file.
    const int CHROM_NOT_FOUND = -1;

    // Process exit status used on fatal errors.
    const int EXIT_FAILURE_STATUS = 1;

    // The very first base of a chromosome (BED coordinates are 0-based).
    const int CHROM_FIRST_BASE = 0;

    // Strand symbol for a feature on the forward strand.
    const char * const STRAND_FORWARD = "+";
    // Strand symbol for a feature on the reverse strand.
    const char * const STRAND_REVERSE = "-";
}


//************************************************
// Class methods and elements
//************************************************
//
// BedFlank implements "bedtools flank": for each feature it creates new
// flanking features to the left and right (and, with -s, honours strand).
// Flank sizes may be given in absolute bp or as a fraction of the feature
// length (-pct). Flanks are clamped to the bounds of the chromosome.
//
class BedFlank {

public:

    // constructor
    BedFlank(string &bedFile, string &genomeFile, bool forceStrand, 
             float leftSlop, float rightSlop, bool fractional,
             bool printHeader);

    // destructor
    ~BedFlank(void);



private:

    string _bedFile;
    string _genomeFile;

    bool   _forceStrand;
    float  _leftFlank;
    float  _rightFlank;
    bool   _fractional;
    bool   _printHeader;

    BedFile *_bed;
    GenomeFile *_genome;

    // methods

    void FlankBed();

    // method to grab requested flank w.r.t. a single BED entry
    void AddFlank(BED &bed, int leftSlop, int rightSlop);
    
    // method to grab requested flank w.r.t. a single BED entry, 
    // while choosing flanks based on strand
    void AddStrandedFlank(BED &bed, int leftSlop, int rightSlop);
};
