/*****************************************************************************
  windowBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef WINDOWBED_H
#define WINDOWBED_H

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
using namespace BamTools;

#include "bedFile.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

//************************************************
// Named constants (replace previously hard-coded
// "magic" literals used throughout windowBed).
//************************************************
namespace WindowBedConstants {
    // Number of BED hit records to pre-allocate per query feature.
    // Chosen as a reasonable upper bound to avoid repeated re-allocation.
    const unsigned int DEFAULT_HITS_RESERVE = 100;

    // When the "A" input is BAM, each alignment is treated as a BED6
    // feature (chrom, start, end, name, score, strand).
    const int BED6_FIELD_COUNT = 6;

    // Minimum overlap fraction passed to the BedFile hit routines.
    // windowBed itself does not enforce a fractional overlap, so it
    // always requests the smallest possible (any) overlap.
    const float ANY_OVERLAP_FRACTION = 0.0;

    // Process exit status used when a fatal error is encountered.
    const int EXIT_FAILURE_STATUS = 1;

    // Strand symbols used when converting BAM alignments to BED.
    const char * const STRAND_FORWARD = "+";
    const char * const STRAND_REVERSE = "-";

    // Read-name suffixes marking the first/second mate of a BAM pair.
    const char * const MATE_ONE_SUFFIX = "/1";
    const char * const MATE_TWO_SUFFIX = "/2";
}

//************************************************
// Class methods and elements
//************************************************
//
// BedWindow implements "bedtools window": for every feature in the
// "A" file it builds a window (by adding slop to the left and right)
// and reports the features in the "B" file that overlap that window.
// The "A" input may be either a BED/GFF/VCF file or a BAM file.
//
class BedWindow {

public:

    // constructor
    BedWindow(string bedAFile, string bedBFile, int leftSlop, int rightSlop,
              bool anyHit, bool noHit, bool writeCount, bool strandWindows,
              bool matchOnSameStrand, bool matchOnDiffStrand, bool bamInput,
              bool bamOutput, bool isUncompressedBam, bool printHeader);

    // destructor
    ~BedWindow(void);

private:

    string _bedAFile;
    string _bedBFile;
    bool _anyHit;
    bool _writeCount;
    int _leftSlop;
    int _rightSlop;
    bool _noHit;
    bool _strandWindows;
    bool _matchOnSameStrand;
    bool _matchOnDiffStrand;
    bool _bamInput;
    bool _bamOutput;
    bool _isUncompressedBam;
    bool _printHeader;

    // instance of a bed file class.
    BedFile *_bedA, *_bedB;

    // methods
    void WindowIntersectBed();
    void WindowIntersectBam(string bamFile);
    void FindWindowOverlaps(const BED &a, vector<BED> &hits);
    bool FindOneOrMoreWindowOverlaps(const BED &a);
    void AddWindow(const BED &a, CHRPOS &fudgeStart, CHRPOS &fudgeEnd);

};
#endif /* WINDOWBED_H */
