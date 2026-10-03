/*****************************************************************************
  pairToBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef INTERSECTBED_H
#define INTERSECTBED_H

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
using namespace BamTools;

#include "bedFile.h"
#include "bedFilePE.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;


//************************************************
// Named constants (replace hard-coded literals).
//************************************************
namespace PairToBedConstants {
    // The nine valid values of the -type option. Using named constants keeps
    // the many string comparisons in the implementation legible and avoids
    // typos drifting apart from the option parser.
    const char * const TYPE_EITHER   = "either";    // overlap on either end
    const char * const TYPE_NEITHER  = "neither";   // no overlap on either end
    const char * const TYPE_BOTH     = "both";      // overlap required on both ends
    const char * const TYPE_XOR      = "xor";       // overlap on exactly one end
    const char * const TYPE_NOTBOTH  = "notboth";   // xor + neither
    const char * const TYPE_ISPAN    = "ispan";     // inner span [end1, start2]
    const char * const TYPE_OSPAN    = "ospan";     // outer span [start1, end2]
    const char * const TYPE_NOTISPAN = "notispan";  // report A if ispan misses B
    const char * const TYPE_NOTOSPAN = "notospan";  // report A if ospan misses B

    // Number of BED hit records to pre-allocate for the plain-text path.
    const unsigned int HITS_RESERVE_BED = 100;
    // Larger reservation used while processing BAM alignment blocks.
    const unsigned int HITS_RESERVE_BAM = 1000;

    // A BAM-derived pair is reported as a full 10-column BEDPE record.
    const int BEDPE_FIELD_COUNT = 10;

    // Minimum overlap fraction requesting "any" (>=1bp) overlap.
    const float ANY_OVERLAP_FRACTION = 0.0;

    // Process exit status used on fatal errors.
    const int EXIT_FAILURE_STATUS = 1;

    // Placeholder chromosome marking an unmapped / absent BEDPE end.
    const char * const UNMAPPED_CHROM = ".";
}



/**************************************************
Helper function protoypes
**************************************************/
void IsCorrectMappingForBEDPE (const BamAlignment &bam, const RefVector &refs, BEDPE &a);



//************************************************
// Class methods and elements
//************************************************
class BedIntersectPE {

public:

    // constructor
    BedIntersectPE(string bedAFilePE, string bedBFile, float overlapFraction,
        string searchType, bool sameStrand, bool diffStrand, bool bamInput, bool bamOutput, bool uncompressedBam, bool useEditDistance);
    // destructor
    ~BedIntersectPE(void);

    void FindOverlaps(const BEDPE &, vector<BED> &hits1, vector<BED> &hits2, const string &type);

    bool FindOneOrMoreOverlaps(const BEDPE &, const string &type);

    void FindSpanningOverlaps(const BEDPE &a, vector<BED> &hits, const string &type);
    bool FindOneOrMoreSpanningOverlaps(const BEDPE &a, const string &type);

    void IntersectBedPE();
    void IntersectBamPE(string bamFile);

    void DetermineBedPEInput();

private:

    string _bedAFilePE;
    string _bedBFile;
    float _overlapFraction;
    string _searchType;
    bool _sameStrand;
    bool _diffStrand;
    bool _useEditDistance;
    bool _bamInput;
    bool _bamOutput;
    bool  _isUncompressedBam;

    // instance of a paired-end bed file class.
    BedFilePE *_bedA;

    // instance of a bed file class.
    BedFile *_bedB;

    inline
    void ConvertBamToBedPE(const BamAlignment &bam1, const BamAlignment &bam2, const RefVector &refs, BEDPE &a) {

        // initialize BEDPE variables
        a.start1 = a.start2 = a.end1 = a.end2 = -1;
        a.chrom1 = a.chrom2 = ".";
        a.strand1 = a.strand2 = '.';
        uint32_t editDistance1, editDistance2;
        editDistance1 = editDistance2 = 0;

        // take the qname from end 1.
        a.name = bam1.Name;

        // end 1
        if (bam1.IsMapped()) {
            a.chrom1  = refs.at(bam1.RefID).RefName;
            a.start1  = bam1.Position;
            a.end1    = bam1.GetEndPosition(false, false);
            a.strand1 = "+";
            if (bam1.IsReverseStrand()) a.strand1 = "-";

            // extract the edit distance from the NM tag
            // if possible. otherwise, complain.
            if (_useEditDistance == true) {
                if (bam1.GetTag("NM", editDistance1) == false) {
                    cerr << "The edit distance tag (NM) was not found in the BAM file.  Please disable -ed.  Exiting\n";
                    exit(1);
                }
            }
        }

        // end 2
        if (bam2.IsMapped()) {
            a.chrom2  = refs.at(bam2.RefID).RefName;
            a.start2  = bam2.Position;
            a.end2    = bam2.GetEndPosition(false, false);
            a.strand2 = "+";
            if (bam2.IsReverseStrand()) a.strand2 = "-";

            // extract the edit distance from the NM tag
            // if possible. otherwise, complain.
            if (_useEditDistance == true) {
                if (bam2.GetTag("NM", editDistance2) == false) {
                    cerr << "The edit distance tag (NM) was not found in the BAM file.  Please disable -ed.  Exiting\n";
                    exit(1);
                }
            }
        }

        // swap the ends if necessary
        if ( a.chrom1 > a.chrom2 || ((a.chrom1 == a.chrom2) && (a.start1 > a.start2)) ) {
            swap(a.chrom1, a.chrom2);
            swap(a.start1, a.start2);
            swap(a.end1, a.end2);
            swap(a.strand1, a.strand2);
        }

        // compute the minimum mapping quality b/w the two ends of the pair.
        a.score = "0";
        if (_useEditDistance == false) {
            if (bam1.IsMapped() == true && bam2.IsMapped() == true)
                a.score = ToString(min(bam1.MapQuality, bam2.MapQuality));
        }
        // BEDPE using edit distance
        else {
            if (bam1.IsMapped() == true && bam2.IsMapped() == true)
                a.score = ToString((int) (editDistance1 + editDistance2));
            else if (bam1.IsMapped() == true)
                a.score = ToString((int) editDistance1);
            else if (bam2.IsMapped() == true)
                a.score = ToString((int) editDistance2);
        }
    };

    inline
    void ProcessBamBlock (const BamAlignment &bam1, const BamAlignment &bam2,
                                          const RefVector &refs,
                                          BamWriter &writer);
};

#endif /* PEINTERSECTBED_H */
