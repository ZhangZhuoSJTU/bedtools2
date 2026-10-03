/*****************************************************************************
  windowBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "windowBed.h"

// Bring the named replacements for previously hard-coded literals into scope.
using namespace WindowBedConstants;


/*
    Constructor

    Stores the user's options, opens the A and B BedFile objects, and
    immediately dispatches to the BED or BAM intersection routine
    depending on whether the A input is BAM.
*/
BedWindow::BedWindow(string bedAFile, string bedBFile, int leftSlop, int rightSlop,
                     bool anyHit, bool noHit, bool writeCount, bool strandWindows,
                     bool matchOnSameStrand, bool matchOnDiffStrand, bool bamInput, 
                     bool bamOutput, bool isUncompressedBam, bool printHeader) {

    _bedAFile      = bedAFile;
    _bedBFile      = bedBFile;

    _leftSlop      = leftSlop;
    _rightSlop     = rightSlop;

    _anyHit              = anyHit;
    _noHit               = noHit;
    _writeCount          = writeCount;
    _strandWindows       = strandWindows;
    _matchOnSameStrand   = matchOnSameStrand;
    _matchOnDiffStrand   = matchOnDiffStrand;
    _bamInput            = bamInput;
    _bamOutput           = bamOutput;
    _isUncompressedBam   = isUncompressedBam;
    _printHeader         = printHeader;

    _bedA          = new BedFile(bedAFile);
    _bedB          = new BedFile(bedBFile);

    if (_bamInput == false)
        WindowIntersectBed();
    else
        WindowIntersectBam(_bedAFile);
}



/*
    Destructor
*/
BedWindow::~BedWindow(void) {
}



void BedWindow::FindWindowOverlaps(const BED &a, vector<BED> &hits) {

    /*
        Adjust the start and end of a based on the requested window
    */

    // update the current feature's start and end
    // according to the slop requested (slop = 0 by default)
    CHRPOS aFudgeStart = 0;
    CHRPOS aFudgeEnd;
    AddWindow(a, aFudgeStart, aFudgeEnd);

    /*
        Now report the hits (if any) based on the window around a.
    */
    // get the hits in B for the A feature
    _bedB->allHits(a.chrom, aFudgeStart, aFudgeEnd, a.strand, hits,
                   _matchOnSameStrand, _matchOnDiffStrand, ANY_OVERLAP_FRACTION, false);

    int numOverlaps = 0;

    // loop through the hits and report those that meet the user's criteria
    vector<BED>::const_iterator h = hits.begin();
    vector<BED>::const_iterator hitsEnd = hits.end();
    for (; h != hitsEnd; ++h) {

        const int s = max(aFudgeStart, h->start);
        const int e = min(aFudgeEnd, h->end);
        const int overlapBases = (e - s);       // the number of overlapping bases b/w a and b
        const int aLength = (a.end - a.start);   // the length of a in b.p.

        if (s < e) {
            // is there enough overlap (default ~ 1bp)
            if ( ((float) overlapBases / (float) aLength) > 0 ) {
                numOverlaps++;
                if (_anyHit == false && _noHit == false && _writeCount == false) {
                    _bedA->reportBedTab(a);
                    _bedB->reportBedNewLine(*h);
                }
            }
        }
    }
    if (_anyHit == true && (numOverlaps >= 1)) {
        _bedA->reportBedNewLine(a); }
    else if (_writeCount == true) {
        _bedA->reportBedTab(a); printf("%d\n", numOverlaps);
    }
    else if (_noHit == true && (numOverlaps == 0)) {
        _bedA->reportBedNewLine(a);
    }
}


bool BedWindow::FindOneOrMoreWindowOverlaps(const BED &a) {

    // update the current feature's start and end
    // according to the slop requested (slop = 0 by default)
    CHRPOS aFudgeStart = 0;
    CHRPOS aFudgeEnd;
    AddWindow(a, aFudgeStart, aFudgeEnd);

    bool overlapsFound = _bedB->anyHits(a.chrom, aFudgeStart, aFudgeEnd, a.strand,
                                        _matchOnSameStrand, _matchOnDiffStrand,
                                        ANY_OVERLAP_FRACTION, false);
    return overlapsFound;
}


// Stream the plain-text "A" BED file, comparing every valid feature
// against the pre-loaded "B" interval map.
void BedWindow::WindowIntersectBed() {

    // load the "B" bed file into a map so
    // that we can easily compare "A" to it for overlaps
    _bedB->loadBedFileIntoMap();

    BED a;
    vector<BED> hits;
    hits.reserve(DEFAULT_HITS_RESERVE);

    _bedA->Open();
    // report A's header first if asked.
    if (_printHeader == true) {
        _bedA->PrintHeader();
    }
    while (_bedA->GetNextBed(a)) {
        if (_bedA->_status == BED_VALID) {
            FindWindowOverlaps(a, hits);
            hits.clear();
        }
    }
    _bedA->Close();
}


// Stream a BAM "A" file: convert every mapped alignment into a BED6
// feature, test it against the "B" map, and emit either BAM or BED.
void BedWindow::WindowIntersectBam(string bamFile) {

    // load the "B" bed file into a map so
    // that we can easily compare "A" to it for overlaps
    _bedB->loadBedFileIntoMap();

    // open the BAM file
    BamReader reader;
    BamWriter writer;
    if (!reader.Open(bamFile)) {
        cerr << "Failed to open BAM file " << bamFile << endl;
        exit(EXIT_FAILURE_STATUS);
    }


    // get header & reference information
    string bamHeader  = reader.GetHeaderText();
    RefVector refs    = reader.GetReferenceData();

    // open a BAM output to stdout if we are writing BAM
    if (_bamOutput == true) {
        // set compression mode
        BamWriter::CompressionMode compressionMode = BamWriter::Compressed;
        if ( _isUncompressedBam ) compressionMode = BamWriter::Uncompressed;
        writer.SetCompressionMode(compressionMode);
        // open our BAM writer
        writer.Open("stdout", bamHeader, refs);
    }

    vector<BED> hits;                   // vector of potential hits
    // reserve some space
    hits.reserve(DEFAULT_HITS_RESERVE);

    _bedA->bedType = BED6_FIELD_COUNT;
    BamAlignment bam;
    bool overlapsFound;
    // get each set of alignments for each pair.
    while (reader.GetNextAlignment(bam)) {

        if (bam.IsMapped()) {
            BED a;
            a.chrom = refs.at(bam.RefID).RefName;
            a.start = bam.Position;
            a.end   = bam.GetEndPosition(false, false);

            // build the name field from the BAM alignment.
            a.name = bam.Name;
            if (bam.IsFirstMate()) a.name += MATE_ONE_SUFFIX;
            if (bam.IsSecondMate()) a.name += MATE_TWO_SUFFIX;

            a.score  = ToString(bam.MapQuality);
            a.strand = STRAND_FORWARD; if (bam.IsReverseStrand()) a.strand = STRAND_REVERSE;

            if (_bamOutput == true) {
                overlapsFound = FindOneOrMoreWindowOverlaps(a);
                if (overlapsFound == true) {
                    if (_noHit == false)
                        writer.SaveAlignment(bam);
                }
                else {
                    if (_noHit == true)
                        writer.SaveAlignment(bam);
                }
            }
            else {
                FindWindowOverlaps(a, hits);
                hits.clear();
            }
        }
        // BAM IsMapped() is false
        else if (_noHit == true) {
            writer.SaveAlignment(bam);
        }
    }

    // close the relevant BAM files.
    reader.Close();
    if (_bamOutput == true) {
        writer.Close();
    }
}


// Compute the windowed [fudgeStart, fudgeEnd) coordinates for feature "a"
// by applying the left/right slop. The start is clamped at 0 so the window
// never extends before the beginning of the chromosome.
void BedWindow::AddWindow(const BED &a, CHRPOS &fudgeStart, CHRPOS &fudgeEnd) {
    // Does the user want to treat the windows based on strand?
    // If so,
    // if "+", then left is left and right is right
    // if "-", the left is right and right is left.
    if (_strandWindows) {
        if (a.strand == "+") {
            if ((int) (a.start - _leftSlop) > 0) 
                fudgeStart = a.start - _leftSlop;
            else fudgeStart = 0;
            fudgeEnd = a.end + _rightSlop;
        }
        else {
            if ((int) (a.start - _rightSlop) > 0) 
                fudgeStart = a.start - _rightSlop;
            else fudgeStart = 0;
            fudgeEnd = a.end + _leftSlop;
        }
    }
    // If not, add the windows irrespective of strand
    else {
        if ((int) (a.start - _leftSlop) > 0) 
            fudgeStart = a.start - _leftSlop;
        else fudgeStart = 0;
        fudgeEnd = a.end + _rightSlop;
    }
}

