/*****************************************************************************
  shuffleBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "bedFile.h"
#include "bedFilePE.h"
#include "GenomeFile.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>
#include <unistd.h>
#include <sys/types.h>
#include <algorithm>  // for binary search
using namespace std;

//************************************************
// Class methods and elements
//************************************************
// BedShuffle randomly repositions the intervals of a BED (or BEDPE) file
// across a genome while preserving each interval's length. The placement
// can be constrained in several ways:
//   * -excl : avoid overlapping a set of excluded regions,
//   * -incl : only land within a set of included regions (size-weighted),
//   * -chrom: keep each interval on its original chromosome,
//   * -noOverlapping: forbid shuffled intervals from overlapping each other.
// A bounded retry loop (-maxTries) guards the constrained modes so that an
// interval that cannot be satisfactorily placed is skipped with a warning.
class BedShuffle {

public:

    // constructor
    BedShuffle(string &bedFile, string &genomeFile, 
               string &excludeFile, string &includeFile, 
               bool haveSeed, bool haveExclude, 
               bool haveInclude, bool sameChrom, 
               float overlapFraction, int seed, 
               bool chooseChrom, bool isBedpe,
               size_t _maxTries, bool noOverlapping,
               bool preventExceedingChromEnd);

    // destructor
    ~BedShuffle(void);

private:

    string _bedFile;
    string _genomeFile;
    string _excludeFile;
    string _includeFile;
    float  _overlapFraction;
    int _seed;
    bool _sameChrom;
    bool _haveExclude;
    bool _haveInclude;
    bool _haveSeed;
    bool _chooseChrom;
    bool _isBedpe;
    size_t _tries;
    size_t _maxTries;
    bool _noOverlapping;
    bool _preventExceedingChromEnd;

    // The BED file from which to compute coverage.
    BedFile *_bed;
    BedFilePE *_bedpe;
    BedFile *_exclude;
    BedFile *_include;

    GenomeFile *_genome;

    vector<string> _chroms;
    int _numChroms;
    vector<string> _includeChroms;
    //    int _numIncludeChroms;
    CHRPOS _genomeSize;

    // include length sum
    long double _cumLen;
    
    // methods
    // Unconstrained shuffle (no -excl / -incl).
    void Shuffle();
    // Shuffle while avoiding the excluded regions (-excl / -noOverlapping).
    void ShuffleWithExclusions();
    // Shuffle while staying inside the included regions (-incl).
    void ShuffleWithInclusions();
    // Shuffle inside included regions while avoiding excluded ones.
    void ShuffleWithInclusionsAndExclusions();

    // Pick a random placement for a single BED interval.
    void ChooseLocus(BED &);
    // Pick a random placement inside a size-weighted -incl interval.
    void ChooseLocusFromInclusionFile(BED &);

    // Pick a random placement for a BEDPE (paired) record.
    void ChoosePairedLocus(BEDPE &b);

    // Emit the standard "could not place after _maxTries" warning. Factored
    // out because the identical message is produced from several shuffle
    // modes when an interval cannot be placed within the retry budget.
    void reportMaxTriesExceeded();
};
