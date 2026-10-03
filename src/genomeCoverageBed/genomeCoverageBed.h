/*****************************************************************************
genomeCoverage.h

(c) 2009 - Aaron Quinlan
Hall Laboratory
Department of Biochemistry and Molecular Genetics
University of Virginia
aaronquinlan@gmail.com

Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "bedFile.h"
#include "NewGenomeFile.h"

#include "BlockedIntervals.h"
#include "api/BamReader.h"
#include "api/BamAux.h"
using namespace BamTools;

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
using namespace std;


//***********************************************
// Named constants (replacing magic literals)
//***********************************************
// genomeCoverageBed reports, per position, how many input intervals
// cover that base of the genome. These constants give names to the
// literals that control how that coverage is emitted.

// BED intervals are 0-based, half-open. The per-base (-d / -dz) report
// can optionally be emitted using 1-based coordinates instead; these
// offsets are added to the internal 0-based position before printing.
static const CHRPOS GENOME_COVERAGE_ZERO_BASED_OFFSET = 0;
static const CHRPOS GENOME_COVERAGE_ONE_BASED_OFFSET  = 1;

// When -scale is left at its default (1.0) the reported depths are
// integral, so output precision/formatting is specialized for this case.
static const float GENOME_COVERAGE_DEFAULT_SCALE = 1.0;

// Sentinels used while walking a chromosome in bedGraph mode to detect
// the first position and transitions between runs of equal depth.
static const CHRPOS GENOME_COVERAGE_NO_START = -1;
static const int    GENOME_COVERAGE_NO_DEPTH = -1;

// The UCSC Genome Browser track definition line emitted for bedGraph
// output when the user also requests a track line (-bg[a] -trackline).
static const char* const GENOME_COVERAGE_BEDGRAPH_TRACK_TYPE =
    "track type=bedGraph";

//***********************************************
// Typedefs
//***********************************************
// depthMap:      depth bin -> number of bases at that depth (per chrom)
// chromDepthMap: chrom name -> that chrom's running coverage depths
typedef map<int, DEPTH, less<int> > depthMap;
typedef map<string, depthMap, less<string> > chromDepthMap;

// histMap:      depth bin -> number of bases at that depth
// chromHistMap: chrom name -> that chrom's depth histogram
typedef map<int, unsigned int, less<int> > histMap;
typedef map<string, histMap, less<string> > chromHistMap;

//************************************************
// Class methods and elements
//************************************************
// BedGenomeCoverage computes the depth of feature coverage across a
// genome from either a BED/GFF/VCF file or a BAM file. Depending on the
// options supplied it can emit one of several report styles:
//   * a genome-wide depth histogram (the default),
//   * a per-base depth report (-d / -dz),
//   * a bedGraph of contiguous equal-depth runs (-bg / -bga).
// Coverage is accumulated one chromosome at a time; the input is assumed
// to be grouped by chromosome so that only one chromosome's worth of
// coverage is held in memory at any moment.
class BedGenomeCoverage {

public:

    // constructor
    BedGenomeCoverage(string bedFile, string genomeFile,
                      bool eachBase, bool startSites,
                      bool bedGraph, bool bedGraphAll,
                      int max, float scale,
                      bool bamInput, bool obeySplits, bool ignoreD,
                      bool filterByStrand, string requestedStrand,
                      bool only_5p_end, bool only_3p_end,
                      bool pair_chip,bool haveSize, int fragmentSize, bool dUTP,
                      bool eachBaseZeroBased,
                      bool add_gb_track_line, string gb_track_line_opts);

    // destructor
    ~BedGenomeCoverage(void);

private:

    // data (parms)
    string _bedFile;
    string _genomeFile;
    bool _bamInput;
    bool _eachBase;
    bool _eachBaseZeroBased;
    bool _startSites;
    bool _bedGraph;
    bool _bedGraphAll;
    int _max;
    float _scale;
    bool _obeySplits;
    bool _ignoreD;
    bool _filterByStrand;
    bool _only_5p_end;
    bool _only_3p_end;
    bool _pair_chip_;
    bool _haveSize;
    bool _dUTP;
    int _fragmentSize;
    bool _add_gb_track_line;
    string _gb_track_line_opts;
    string _requestedStrand;

    BedFile *_bed;
    NewGenomeFile *_genome;

    // data for internal processing
    chromDepthMap _chromCov;
    string _currChromName ;
    vector<DEPTH> _currChromCoverage;
    chromHistMap _currChromDepthHist;
    CHRPOS _currChromSize ;
    set<string> _visitedChromosomes;


    // methods
    // Drive coverage calculation from a BED/GFF/VCF or a BAM input.
    void CoverageBed();
    void CoverageBam(string bamFile);
    void LoadBamHeaderIntoGenomeFile(const string &bamFile);
    // Emit the coverage for a single chromosome in the requested style.
    void ReportChromCoverage(const vector<DEPTH> &, const CHRPOS &chromSize, const string &chrom, chromHistMap&);
    // Emit the aggregate depth histogram across the whole genome.
    void ReportGenomeCoverage(chromHistMap &chromDepthHist);
    // Emit a chromosome's coverage as contiguous equal-depth bedGraph runs.
    void ReportChromCoverageBedGraph(const vector<DEPTH> &chromCov, const CHRPOS &chromSize, const string &chrom);
    // Clear the per-chromosome accumulator between chromosomes.
    void ResetChromCoverage();
    // Begin accumulating coverage for a newly encountered chromosome.
    void StartNewChrom (const string& chrom);
    // Record that an interval spans [start, end] on the current chromosome.
    void AddCoverage (CHRPOS start, CHRPOS end);
    // Record coverage for each block of a split/blocked feature.
    void AddBlockedCoverage(const vector<BED> &bedBlocks);
    // Emit the genome-wide histogram once all chromosomes are processed.
    void PrintFinalCoverage();
    // Emit zero-coverage rows for chromosomes absent from the input.
    void PrintEmptyChromosomes();
    // Emit the optional UCSC browser track line for bedGraph output.
    void PrintTrackDefinitionLine();
};
