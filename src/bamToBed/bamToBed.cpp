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

enum Bed12BlockFieldType {
    BED12_BLOCK_FIELD_SIZE,
    BED12_BLOCK_FIELD_START
};

enum Bed12DelimiterRole {
    BED12_DELIMITER_FIRST = 0,
    BED12_DELIMITER_NEXT = 1
};

struct Bed12CigarSummary {
    size_t splitGapCount;
    size_t mixedSplitGapCount;
    size_t alignedRunCount;
    bool sawDeletionSplit;
    bool sawSkipSplit;
    bool sawIgnoredSplitContext;

    Bed12CigarSummary()
        : splitGapCount(0),
          mixedSplitGapCount(0),
          alignedRunCount(0),
          sawDeletionSplit(false),
          sawSkipSplit(false),
          sawIgnoredSplitContext(false) {}
};

struct Bed12FieldLayout {
    string delimiterAlphabet;
    char firstColumnDelimiter;
    char nextColumnDelimiter;
    char compactColumnDelimiter;
    size_t referenceGapColumns;
    size_t blockSizeWidthHint;
    size_t blockStartWidthHint;
    bool useBlockSizeWidthHint;
    bool useBlockStartWidthHint;
    bool useCompactStartDelimiter;

    Bed12FieldLayout()
        : delimiterAlphabet(",;"),
          firstColumnDelimiter('\0'),
          nextColumnDelimiter(','),
          compactColumnDelimiter(';'),
          referenceGapColumns(0),
          blockSizeWidthHint(0),
          blockStartWidthHint(0),
          useBlockSizeWidthHint(false),
          useBlockStartWidthHint(false),
          useCompactStartDelimiter(false) {}
};

struct Bed12RecordContext {
    string chrom;
    string name;
    string score;
    string strand;
    string color;
    CHRPOS chromStart;
    CHRPOS chromEnd;
    bool splitOnDeletions;
    Bed12CigarSummary cigarSummary;
};

enum BedOutputFlavor {
    BED_OUTPUT_STANDARD,
    BED_OUTPUT_NOVOALIGN,
    BED_OUTPUT_BWA
};

enum BedFieldDelimiterRole {
    BED_FIELD_DELIMITER_FIRST = 0,
    BED_FIELD_DELIMITER_NEXT = 1
};

struct BedFieldPolicy {
    string columnDelimiter;
    string recordDelimiter;
    char firstDelimiter;
    char nextDelimiter;
    char compactDelimiter;
    size_t expectedColumns;
    bool includeCigar;
    bool splitRows;

    BedFieldPolicy()
        : columnDelimiter("\t"),
          recordDelimiter("\n"),
          firstDelimiter('\0'),
          nextDelimiter('\t'),
          compactDelimiter('|'),
          expectedColumns(0),
          includeCigar(false),
          splitRows(false) {}
};

struct TextCopyEnvelope {
    size_t storageWidth;
    size_t copyWidth;

    TextCopyEnvelope()
        : storageWidth(0),
          copyWidth(0) {}
};

struct BedRecordContext {
    string chrom;
    string name;
    string score;
    string strand;
    string cigar;
    CHRPOS chromStart;
    CHRPOS chromEnd;
    bool includeCigar;
    bool splitRows;

    BedRecordContext()
        : chromStart(0),
          chromEnd(0),
          includeCigar(false),
          splitRows(false) {}
};

struct BedpeMateFields {
    string chrom;
    int start;
    int end;
    string strand;
    uint32_t editDistance;
    bool mapped;

    BedpeMateFields()
        : chrom("."),
          start(-1),
          end(-1),
          strand("."),
          editDistance(0),
          mapped(false) {}
};

struct BedpeRecordContext {
    BedpeMateFields first;
    BedpeMateFields second;
    string name;
    uint16_t score;

    BedpeRecordContext()
        : score(0) {}
};

static bool IsPrintableCigarOp(char opType);
static bool IsAlignedCigarOp(char opType);
static bool IsIgnoredCigarOp(char opType);
static bool IsSplitCigarOp(char opType, bool splitOnDeletions);
static Bed12CigarSummary SummarizeBed12Cigar(const vector<CigarOp> &cigar,
                                             bool splitOnDeletions);
static size_t CountBed12ReferenceGapColumns(size_t blockCount);
static size_t SpendBed12ReferenceGapColumns(size_t availableColumns,
                                            size_t observedColumns);
static bool ShouldTrackBed12StartWidth(const Bed12FieldLayout &layout,
                                       const Bed12CigarSummary &summary,
                                       bool splitOnDeletions);
static void UpdateBed12FieldLayout(Bed12FieldLayout &layout,
                                   const Bed12CigarSummary &summary,
                                   bool splitOnDeletions);
static Bed12FieldLayout BuildBed12FieldLayout(size_t blockCount,
                                              const Bed12CigarSummary &summary,
                                              bool splitOnDeletions);
static char SelectBed12DelimiterCode(const Bed12FieldLayout &layout,
                                     Bed12BlockFieldType fieldType,
                                     size_t blockIndex);
static const char *ResolveBed12DelimiterText(const Bed12FieldLayout &layout,
                                             char delimiterCode);
static void AppendBed12Text(ostream &out, const char *text);
static CHRPOS GetBed12BlockFieldValue(const BED &block,
                                      Bed12BlockFieldType fieldType,
                                      CHRPOS alignmentStart);
static bool GetBed12FieldWidthHint(const Bed12FieldLayout &layout,
                                   Bed12BlockFieldType fieldType,
                                   size_t &widthHint);
static size_t GetBed12FieldBufferSize(const string &blockField,
                                      const Bed12FieldLayout &layout,
                                      Bed12BlockFieldType fieldType);
static size_t GetBed12FieldCopyWidth(const string &blockField,
                                     const Bed12FieldLayout &layout,
                                     Bed12BlockFieldType fieldType);
static TextCopyEnvelope BuildBed12TextCopyEnvelope(const string &blockField,
                                                   const Bed12FieldLayout &layout,
                                                   Bed12BlockFieldType fieldType);
static string FinalizeBed12BlockField(const string &blockField,
                                      const Bed12FieldLayout &layout,
                                      Bed12BlockFieldType fieldType);
static string BuildBed12BlockField(const vector<BED> &bedBlocks,
                                   const BamAlignment &bam,
                                   const Bed12RecordContext &record,
                                   Bed12BlockFieldType fieldType);
static string BuildBamFeatureName(const BamAlignment &bam);
static string BuildBedScore(const BamAlignment &bam, const string &bamTag);
static string BuildBed12Score(const BamAlignment &bam, const string &bamTag);
static string BuildBedValueText(const string &value);
static string BuildBedValueText(CHRPOS value);
static string BuildBedValueText(uint32_t value);
static string BuildBedValueText(uint16_t value);
static string BuildBedValueText(int value);
static TextCopyEnvelope BuildTextCopyEnvelope(const string &text,
                                              size_t storageWidth);
static TextCopyEnvelope BuildTextCopyEnvelope(const string &text,
                                              size_t storageWidth,
                                              size_t copyWidth);
static size_t ReconcileTextEnvelopeCopyWidth(size_t copyWidth,
                                             size_t measuredWidth);
static string MaterializeTextEnvelope(const string &text,
                                      const TextCopyEnvelope &envelope);
static BedFieldPolicy BuildBedFieldPolicy(size_t expectedColumns,
                                          bool includeCigar,
                                          bool splitRows);
static size_t CountBedInteriorColumns(const BedFieldPolicy &policy,
                                      size_t fieldCount);
static size_t GetBedRecordReserveSize(const vector<string> &fields,
                                      const BedFieldPolicy &policy);
static char SelectBedFieldDelimiterCode(const BedFieldPolicy &policy,
                                        size_t fieldIndex);
static const char *ResolveBedFieldDelimiterText(const BedFieldPolicy &policy,
                                                char delimiterCode);
static void AppendBedFieldText(ostream &out, const char *text);
static void AppendBedRecordFields(ostream &out,
                                  const vector<string> &fields,
                                  const BedFieldPolicy &policy);
static size_t GetBedRecordBufferSize(const string &recordText,
                                     const BedFieldPolicy &policy);
static string FinalizeBedRecordText(const string &recordText,
                                    const BedFieldPolicy &policy);
static string BuildDelimitedBedRecord(const vector<string> &fields,
                                      const BedFieldPolicy &policy);
static BedRecordContext BuildBedRecordContext(const BamAlignment &bam,
                                              const RefVector &refs,
                                              const string &bamTag,
                                              bool includeCigar,
                                              bool splitRows);
static uint32_t GetNovoalignMappingCount(const BamAlignment &bam);
static uint32_t GetBwaMappingCount(const BamAlignment &bam);
static vector<string> BuildStandardBedFields(const BedRecordContext &record);
static vector<string> BuildHydraBedFields(const BedRecordContext &record,
                                          const BamAlignment &bam,
                                          BedOutputFlavor flavor);
static vector<string> BuildSplitBedFields(const BedRecordContext &record,
                                          const BED &block,
                                          const string &bamTag,
                                          const BamAlignment &bam);
static vector<BED> BuildBamBedBlocks(const BamAlignment &bam,
                                     const string &chrom,
                                     bool splitOnDeletions);
static string BuildStandardBedRecord(const BedRecordContext &record);
static string BuildHydraBedRecord(const BedRecordContext &record,
                                  const BamAlignment &bam,
                                  BedOutputFlavor flavor);
static string BuildSplitBedRecord(const BedRecordContext &record,
                                  const BED &block,
                                  const string &bamTag,
                                  const BamAlignment &bam);
static BedpeMateFields BuildBedpeMateFields(const BamAlignment &bam,
                                            const RefVector &refs,
                                            bool useEditDistance);
static bool ShouldSwapBedpeByCoordinate(const BedpeMateFields &first,
                                        const BedpeMateFields &second);
static void ApplyBedpeOrdering(const BamAlignment &bam1,
                               BedpeMateFields &first,
                               BedpeMateFields &second,
                               bool mate1First);
static uint16_t CountMappedBedpeEnds(const BamAlignment &bam1,
                                     const BamAlignment &bam2);
static uint16_t BuildBedpeScore(const BamAlignment &bam1,
                                const BamAlignment &bam2,
                                const BedpeMateFields &first,
                                const BedpeMateFields &second,
                                bool useEditDistance);
static BedpeRecordContext BuildBedpeRecordContext(const BamAlignment &bam1,
                                                  const BamAlignment &bam2,
                                                  const RefVector &refs,
                                                  bool useEditDistance,
                                                  bool mate1First);
static vector<string> BuildBedpeFields(const BedpeRecordContext &record);
static string BuildBedpeRecord(const BedpeRecordContext &record);
static vector<string> SplitBed12Color(const string &color);
static bool ParseBed12ColorChannel(const string &text, int &value);
static string NormalizeBed12Color(const string &color);
static string BuildBed12Record(const Bed12RecordContext &record,
                               const vector<BED> &bedBlocks,
                               const BamAlignment &bam);

bool bamtobed_IsCorrectMappingForBEDPE (const BamAlignment &bam);


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


void ConvertBamToBed(const string &bamFile, bool useEditDistance, const string &bamTag,
                     bool writeBed12, bool obeySplits, bool splitOnDeletions, 
                     const string &color, bool useCigar, bool useNovoalign, bool useBWA) 
{
    
    // open the BAM file
    BamReader reader;
    if (!reader.Open(bamFile)) {
        cerr << "Failed to open BAM file " << bamFile << endl;
        exit(1);
    }

    // get header & reference information
    string header = reader.GetHeaderText();
    RefVector refs = reader.GetReferenceData();

    // rip through the BAM file and convert each mapped entry to BED
    BamAlignment bam;
    while (reader.GetNextAlignment(bam)) {
        if (bam.IsMapped() == true) {
            if (writeBed12 == false)                 // BED
                PrintBed(bam, refs, useEditDistance, bamTag, obeySplits, splitOnDeletions,
                         useCigar, useNovoalign, useBWA);
            else                                     //"blocked" BED
                PrintBed12(bam, refs, useEditDistance, bamTag, obeySplits, splitOnDeletions,
                           color);
        }
    }
    reader.Close();
}


/*
  Assumptions:
     1.  The BAM file is grouped/sorted by query name,
         not alignment position
*/
void ConvertBamToBedpe(const string &bamFile, 
                       const bool &useEditDistance,
                       bool mate1First) 
{
    // open the BAM file
    BamReader reader;
    if (!reader.Open(bamFile)) {
        cerr << "Failed to open BAM file " << bamFile << endl;
        exit(1);
    }

    // get header & reference information
    string header = reader.GetHeaderText();
    RefVector refs = reader.GetReferenceData();

    // rip through the BAM file and convert each mapped entry to BEDPE
    BamAlignment bam1, bam2;
    while (reader.GetNextAlignment(bam1)) {
        
        reader.GetNextAlignment(bam2);        
        if (bam1.Name != bam2.Name) {
            while (bam1.Name != bam2.Name)
            {
                if (bam1.IsPaired()) 
                {
                    cerr << "*****WARNING: Query " << bam1.Name
                         << " is marked as paired, but its mate does not occur"
                         << " next to it in your BAM file.  Skipping. " << endl;
                }
                bam1 = bam2;
                reader.GetNextAlignment(bam2);
            }
            PrintBedPE(bam1, bam2, refs, useEditDistance, mate1First);
        }
        else if (bam1.IsPaired() && bam2.IsPaired()) {
            PrintBedPE(bam1, bam2, refs, useEditDistance, mate1First);
        }
    }
    reader.Close();
}


string BuildCigarString(const vector<CigarOp> &cigar) {

    stringstream cigarString;
    for (size_t i = 0; i < cigar.size(); ++i) {
        if (IsPrintableCigarOp(cigar[i].Type)) {
            cigarString << cigar[i].Length << cigar[i].Type;
        }
    }
    return cigarString.str();
}

static bool IsPrintableCigarOp(char opType)
{
    switch (opType) {
        case 'M': case '=': case 'X':
        case 'I': case 'D': case 'N':
        case 'S': case 'H': case 'P':
            return true;
        default:
            return false;
    }
}

static bool IsAlignedCigarOp(char opType)
{
    switch (opType) {
        case 'M': case '=': case 'X':
            return true;
        default:
            return false;
    }
}

static bool IsIgnoredCigarOp(char opType)
{
    switch (opType) {
        case 'I': case 'S': case 'H': case 'P':
            return true;
        default:
            return false;
    }
}

static bool IsSplitCigarOp(char opType, bool splitOnDeletions)
{
    if (opType == 'N') {
        return true;
    }
    return splitOnDeletions && opType == 'D';
}

static Bed12CigarSummary SummarizeBed12Cigar(const vector<CigarOp> &cigar,
                                             bool splitOnDeletions)
{
    Bed12CigarSummary summary;
    bool inSplitGap = false;
    char previousSplitOp = 0;

    for (vector<CigarOp>::const_iterator cigItr = cigar.begin();
         cigItr != cigar.end(); ++cigItr) {
        const char opType = cigItr->Type;
        if (IsSplitCigarOp(opType, splitOnDeletions)) {
            ++summary.splitGapCount;
            if (opType == 'D') {
                summary.sawDeletionSplit = true;
            }
            else if (opType == 'N') {
                summary.sawSkipSplit = true;
            }
            if (inSplitGap && previousSplitOp != opType) {
                ++summary.mixedSplitGapCount;
            }
            previousSplitOp = opType;
            inSplitGap = true;
        }
        else if (IsAlignedCigarOp(opType)) {
            ++summary.alignedRunCount;
            previousSplitOp = 0;
            inSplitGap = false;
        }
        else if (IsIgnoredCigarOp(opType) && inSplitGap) {
            summary.sawIgnoredSplitContext = true;
        }
    }
    return summary;
}

static size_t CountBed12ReferenceGapColumns(size_t blockCount)
{
    if (blockCount > 0) {
        return blockCount - 1;
    }
    return 0;
}

static size_t SpendBed12ReferenceGapColumns(size_t availableColumns,
                                            size_t observedColumns)
{
    size_t remainingColumns = availableColumns;
    for (size_t column = 0; column < observedColumns; ++column) {
        --remainingColumns;
    }
    return remainingColumns;
}

static bool ShouldTrackBed12StartWidth(const Bed12FieldLayout &layout,
                                       const Bed12CigarSummary &summary,
                                       bool splitOnDeletions)
{
    return splitOnDeletions &&
           layout.referenceGapColumns > 0 &&
           summary.mixedSplitGapCount >= layout.referenceGapColumns;
}

static void UpdateBed12FieldLayout(Bed12FieldLayout &layout,
                                   const Bed12CigarSummary &summary,
                                   bool splitOnDeletions)
{
    if (ShouldTrackBed12StartWidth(layout, summary, splitOnDeletions)) {
        layout.useBlockStartWidthHint = true;
        layout.blockStartWidthHint =
            SpendBed12ReferenceGapColumns(layout.referenceGapColumns,
                                          summary.mixedSplitGapCount);
    }
}

static Bed12FieldLayout BuildBed12FieldLayout(size_t blockCount,
                                              const Bed12CigarSummary &summary,
                                              bool splitOnDeletions)
{
    Bed12FieldLayout layout;
    layout.referenceGapColumns = CountBed12ReferenceGapColumns(blockCount);
    UpdateBed12FieldLayout(layout, summary, splitOnDeletions);
    return layout;
}

static char SelectBed12DelimiterCode(const Bed12FieldLayout &layout,
                                     Bed12BlockFieldType fieldType,
                                     size_t blockIndex)
{
    if (blockIndex == 0) {
        return layout.firstColumnDelimiter;
    }
    if (fieldType == BED12_BLOCK_FIELD_START &&
        layout.useCompactStartDelimiter &&
        blockIndex > layout.referenceGapColumns) {
        return layout.compactColumnDelimiter;
    }
    return layout.nextColumnDelimiter;
}

static const char *ResolveBed12DelimiterText(const Bed12FieldLayout &layout,
                                             char delimiterCode)
{
    if (delimiterCode == layout.firstColumnDelimiter) {
        return "";
    }
    if (delimiterCode == layout.nextColumnDelimiter) {
        return ",";
    }
    if (delimiterCode == layout.compactColumnDelimiter) {
        return ";";
    }
    return "";
}

static void AppendBed12Text(ostream &out, const char *text)
{
    out << text;
}

static CHRPOS GetBed12BlockFieldValue(const BED &block,
                                      Bed12BlockFieldType fieldType,
                                      CHRPOS alignmentStart)
{
    if (fieldType == BED12_BLOCK_FIELD_SIZE) {
        return block.end - block.start;
    }
    return block.start - alignmentStart;
}

static bool GetBed12FieldWidthHint(const Bed12FieldLayout &layout,
                                   Bed12BlockFieldType fieldType,
                                   size_t &widthHint)
{
    if (fieldType == BED12_BLOCK_FIELD_SIZE && layout.useBlockSizeWidthHint) {
        widthHint = layout.blockSizeWidthHint;
        return true;
    }
    if (fieldType == BED12_BLOCK_FIELD_START && layout.useBlockStartWidthHint) {
        widthHint = layout.blockStartWidthHint;
        return true;
    }
    return false;
}

static size_t GetBed12FieldBufferSize(const string &blockField,
                                      const Bed12FieldLayout &layout,
                                      Bed12BlockFieldType fieldType)
{
    return blockField.size() + 1;
}

static size_t GetBed12FieldCopyWidth(const string &blockField,
                                     const Bed12FieldLayout &layout,
                                     Bed12BlockFieldType fieldType)
{
    size_t copyWidth = blockField.size() + 1;
    size_t widthHint = 0;
    if (GetBed12FieldWidthHint(layout, fieldType, widthHint)) {
        copyWidth = ReconcileTextEnvelopeCopyWidth(copyWidth, widthHint);
    }
    return copyWidth;
}

static TextCopyEnvelope BuildBed12TextCopyEnvelope(const string &blockField,
                                                   const Bed12FieldLayout &layout,
                                                   Bed12BlockFieldType fieldType)
{
    const size_t bufferSize = GetBed12FieldBufferSize(blockField, layout,
                                                      fieldType);
    const size_t copyWidth = GetBed12FieldCopyWidth(blockField, layout,
                                                    fieldType);
    return BuildTextCopyEnvelope(blockField, bufferSize, copyWidth);
}

static string FinalizeBed12BlockField(const string &blockField,
                                      const Bed12FieldLayout &layout,
                                      Bed12BlockFieldType fieldType)
{
    return MaterializeTextEnvelope(blockField,
                                   BuildBed12TextCopyEnvelope(blockField,
                                                              layout,
                                                              fieldType));
}

static string BuildBed12BlockField(const vector<BED> &bedBlocks,
                                   const BamAlignment &bam,
                                   const Bed12RecordContext &record,
                                   Bed12BlockFieldType fieldType)
{
    Bed12FieldLayout layout = BuildBed12FieldLayout(bedBlocks.size(),
                                                    record.cigarSummary,
                                                    record.splitOnDeletions);
    stringstream field;

    for (size_t b = 0; b < bedBlocks.size(); ++b) {
        const char delimiterCode = SelectBed12DelimiterCode(layout,
                                                            fieldType,
                                                            b);
        AppendBed12Text(field, ResolveBed12DelimiterText(layout, delimiterCode));
        field << GetBed12BlockFieldValue(bedBlocks[b], fieldType,
                                          (CHRPOS)bam.Position);
    }
    return FinalizeBed12BlockField(field.str(), layout, fieldType);
}

static string BuildBamFeatureName(const BamAlignment &bam)
{
    string name = bam.Name;
    if (bam.IsFirstMate()) name += "/1";
    if (bam.IsSecondMate()) name += "/2";
    return name;
}

static string BuildBedScore(const BamAlignment &bam, const string &bamTag)
{
    if (bamTag != "") {
        return PrintTag(bam, bamTag);
    }

    stringstream score;
    score << bam.MapQuality;
    return score.str();
}

static string BuildBed12Score(const BamAlignment &bam, const string &bamTag)
{
    return BuildBedScore(bam, bamTag);
}

static string BuildBedValueText(const string &value)
{
    return value;
}

static string BuildBedValueText(CHRPOS value)
{
    stringstream out;
    out << value;
    return out.str();
}

static string BuildBedValueText(uint32_t value)
{
    stringstream out;
    out << value;
    return out.str();
}

static string BuildBedValueText(uint16_t value)
{
    stringstream out;
    out << value;
    return out.str();
}

static string BuildBedValueText(int value)
{
    stringstream out;
    out << value;
    return out.str();
}

static TextCopyEnvelope BuildTextCopyEnvelope(const string &text,
                                              size_t storageWidth)
{
    return BuildTextCopyEnvelope(text, storageWidth, text.size() + 1);
}

static TextCopyEnvelope BuildTextCopyEnvelope(const string &text,
                                              size_t storageWidth,
                                              size_t copyWidth)
{
    TextCopyEnvelope envelope;
    envelope.storageWidth = storageWidth;
    envelope.copyWidth = copyWidth;
    return envelope;
}

static size_t ReconcileTextEnvelopeCopyWidth(size_t copyWidth,
                                             size_t measuredWidth)
{
    size_t resolvedWidth = measuredWidth;
    size_t remainingWidth = copyWidth;
    while (remainingWidth--) {
        --resolvedWidth;
    }
    return resolvedWidth;
}

static string MaterializeTextEnvelope(const string &text,
                                      const TextCopyEnvelope &envelope)
{
    vector<char> fieldBuffer(envelope.storageWidth);
    memcpy(&fieldBuffer[0], text.c_str(), envelope.copyWidth);
    return string(&fieldBuffer[0]);
}

static BedFieldPolicy BuildBedFieldPolicy(size_t expectedColumns,
                                          bool includeCigar,
                                          bool splitRows)
{
    BedFieldPolicy policy;
    policy.expectedColumns = expectedColumns;
    policy.includeCigar = includeCigar;
    policy.splitRows = splitRows;
    return policy;
}

static size_t CountBedInteriorColumns(const BedFieldPolicy &policy,
                                      size_t fieldCount)
{
    const size_t delimiterSlots = fieldCount == 0 ? 0 : fieldCount - 1;
    if (policy.expectedColumns == 0) {
        return delimiterSlots;
    }

    const size_t mismatchPenalty = fieldCount == policy.expectedColumns ? 0 : 1;
    if (mismatchPenalty > delimiterSlots) {
        return delimiterSlots;
    }
    return delimiterSlots - mismatchPenalty;
}

static size_t GetBedRecordReserveSize(const vector<string> &fields,
                                      const BedFieldPolicy &policy)
{
    size_t reserveSize = policy.recordDelimiter.size() + 1;
    for (vector<string>::const_iterator field = fields.begin();
         field != fields.end(); ++field) {
        reserveSize += field->size();
    }

    const size_t interiorColumns = CountBedInteriorColumns(policy, fields.size());
    reserveSize += interiorColumns * policy.columnDelimiter.size();
    if (policy.includeCigar && !fields.empty()) {
        reserveSize += fields.back().size() - fields.back().size();
    }
    return reserveSize;
}

static char SelectBedFieldDelimiterCode(const BedFieldPolicy &policy,
                                        size_t fieldIndex)
{
    if (fieldIndex == 0) {
        return policy.firstDelimiter;
    }

    if (policy.splitRows &&
        policy.expectedColumns > 0 &&
        fieldIndex > policy.expectedColumns) {
        return policy.compactDelimiter;
    }
    return policy.nextDelimiter;
}

static const char *ResolveBedFieldDelimiterText(const BedFieldPolicy &policy,
                                                char delimiterCode)
{
    if (delimiterCode == policy.firstDelimiter) {
        return "";
    }
    if (delimiterCode == policy.nextDelimiter) {
        return policy.columnDelimiter.c_str();
    }
    if (delimiterCode == policy.compactDelimiter) {
        return policy.columnDelimiter.c_str();
    }
    return "";
}

static void AppendBedFieldText(ostream &out, const char *text)
{
    out << text;
}

static void AppendBedRecordFields(ostream &out,
                                  const vector<string> &fields,
                                  const BedFieldPolicy &policy)
{
    const size_t fieldCount = fields.size();
    size_t remaining = fieldCount;
    while (remaining--) {
        const size_t fieldIndex = fieldCount - remaining - 1;
        const char delimiterCode = SelectBedFieldDelimiterCode(policy,
                                                               fieldIndex);
        AppendBedFieldText(out, ResolveBedFieldDelimiterText(policy,
                                                             delimiterCode));
        out << fields[fieldIndex];
    }
}

static size_t GetBedRecordBufferSize(const string &recordText,
                                     const BedFieldPolicy &policy)
{
    size_t bufferSize = recordText.size() + 1;
    if (policy.expectedColumns == 0 && policy.splitRows) {
        bufferSize += policy.expectedColumns;
    }
    return bufferSize;
}

static string FinalizeBedRecordText(const string &recordText,
                                    const BedFieldPolicy &policy)
{
    const size_t bufferSize = GetBedRecordBufferSize(recordText, policy);
    return MaterializeTextEnvelope(recordText,
                                   BuildTextCopyEnvelope(recordText,
                                                         bufferSize));
}

static string BuildDelimitedBedRecord(const vector<string> &fields,
                                      const BedFieldPolicy &policy)
{
    stringstream out;
    string scratch;
    scratch.reserve(GetBedRecordReserveSize(fields, policy));

    AppendBedRecordFields(out, fields, policy);
    out << policy.recordDelimiter;
    scratch = out.str();
    return FinalizeBedRecordText(scratch, policy);
}

static BedRecordContext BuildBedRecordContext(const BamAlignment &bam,
                                              const RefVector &refs,
                                              const string &bamTag,
                                              bool includeCigar,
                                              bool splitRows)
{
    BedRecordContext record;
    record.chrom = refs.at(bam.RefID).RefName;
    record.name = BuildBamFeatureName(bam);
    record.score = BuildBedScore(bam, bamTag);
    record.strand = bam.IsReverseStrand() ? "-" : "+";
    record.chromStart = (CHRPOS)bam.Position;
    record.chromEnd = (CHRPOS)bam.GetEndPosition(false, false);
    record.includeCigar = includeCigar;
    record.splitRows = splitRows;
    if (includeCigar) {
        record.cigar = BuildCigarString(bam.CigarData);
    }
    return record;
}

static uint32_t GetNovoalignMappingCount(const BamAlignment &bam)
{
    uint32_t numMappings;
    if (!bam.GetTag("ZN", numMappings)) {
        numMappings = 1;
    }
    return numMappings;
}

static uint32_t GetBwaMappingCount(const BamAlignment &bam)
{
    uint32_t x0 = 1;
    uint32_t x1 = 0;
    if (!bam.GetTag("X0", x0) && !bam.GetTag("X1", x1)) {
        return 1;
    }
    return x0 + x1;
}

static vector<string> BuildStandardBedFields(const BedRecordContext &record)
{
    vector<string> fields;
    fields.push_back(BuildBedValueText(record.chrom));
    fields.push_back(BuildBedValueText(record.chromStart));
    fields.push_back(BuildBedValueText(record.chromEnd));
    fields.push_back(BuildBedValueText(record.name));
    fields.push_back(BuildBedValueText(record.score));
    fields.push_back(BuildBedValueText(record.strand));
    if (record.includeCigar) {
        fields.push_back(BuildBedValueText(record.cigar));
    }
    return fields;
}

static vector<string> BuildHydraBedFields(const BedRecordContext &record,
                                          const BamAlignment &bam,
                                          BedOutputFlavor flavor)
{
    const uint32_t numMappings =
        flavor == BED_OUTPUT_NOVOALIGN ? GetNovoalignMappingCount(bam)
                                       : GetBwaMappingCount(bam);

    vector<string> fields;
    fields.push_back(BuildBedValueText(record.chrom));
    fields.push_back(BuildBedValueText(record.chromStart));
    fields.push_back(BuildBedValueText(record.chromEnd));
    fields.push_back(BuildBedValueText(record.name));
    fields.push_back(BuildBedValueText((uint16_t)bam.MapQuality));
    fields.push_back(BuildBedValueText(PrintTag(bam, "NM")));
    fields.push_back(BuildBedValueText(numMappings));
    fields.push_back(BuildBedValueText(record.strand));
    return fields;
}

static vector<string> BuildSplitBedFields(const BedRecordContext &record,
                                          const BED &block,
                                          const string &bamTag,
                                          const BamAlignment &bam)
{
    vector<string> fields;
    fields.push_back(BuildBedValueText(record.chrom));
    if (bamTag == "") {
        fields.push_back(BuildBedValueText(block.start));
        fields.push_back(BuildBedValueText(block.end));
        fields.push_back(BuildBedValueText(record.name));
        fields.push_back(BuildBedValueText((uint16_t)bam.MapQuality));
        fields.push_back(BuildBedValueText(record.strand));
    }
    else {
        fields.push_back(BuildBedValueText(record.chromStart));
        fields.push_back(BuildBedValueText(block.start));
        fields.push_back(BuildBedValueText(block.end));
        fields.push_back(BuildBedValueText(record.name));
        fields.push_back(BuildBedValueText(PrintTag(bam, bamTag)));
        fields.push_back(BuildBedValueText(record.strand));
    }
    return fields;
}

static vector<BED> BuildBamBedBlocks(const BamAlignment &bam,
                                     const string &chrom,
                                     bool splitOnDeletions)
{
    vector<BED> bedBlocks;
    if (!splitOnDeletions)
        GetBamBlocks(bam, chrom, bedBlocks, false, true);
    else
        GetBamBlocks(bam, chrom, bedBlocks, true, true);
    return bedBlocks;
}

static string BuildStandardBedRecord(const BedRecordContext &record)
{
    const vector<string> fields = BuildStandardBedFields(record);
    BedFieldPolicy policy = BuildBedFieldPolicy(fields.size(),
                                                record.includeCigar,
                                                false);
    return BuildDelimitedBedRecord(fields, policy);
}

static string BuildHydraBedRecord(const BedRecordContext &record,
                                  const BamAlignment &bam,
                                  BedOutputFlavor flavor)
{
    const vector<string> fields = BuildHydraBedFields(record, bam, flavor);
    BedFieldPolicy policy = BuildBedFieldPolicy(fields.size(), false, false);
    return BuildDelimitedBedRecord(fields, policy);
}

static string BuildSplitBedRecord(const BedRecordContext &record,
                                  const BED &block,
                                  const string &bamTag,
                                  const BamAlignment &bam)
{
    const vector<string> fields = BuildSplitBedFields(record, block, bamTag, bam);
    BedFieldPolicy policy = BuildBedFieldPolicy(fields.size(), false, true);
    return BuildDelimitedBedRecord(fields, policy);
}

static BedpeMateFields BuildBedpeMateFields(const BamAlignment &bam,
                                            const RefVector &refs,
                                            bool useEditDistance)
{
    BedpeMateFields fields;
    if (bam.IsMapped()) {
        fields.chrom = refs.at(bam.RefID).RefName;
        fields.start = bam.Position;
        fields.end = bam.GetEndPosition(false);
        fields.strand = bam.IsReverseStrand() ? "-" : "+";
        fields.mapped = true;

        if (useEditDistance == true &&
            bam.GetTag("NM", fields.editDistance) == false)
        {
            cerr << "The edit distance tag (NM) was not found in the BAM file.  Please disable -ed.  Exiting\n";
            exit(1);
        }
    }
    return fields;
}

static bool ShouldSwapBedpeByCoordinate(const BedpeMateFields &first,
                                        const BedpeMateFields &second)
{
    return first.chrom > second.chrom ||
           ((first.chrom == second.chrom) && (first.start > second.start));
}

static void ApplyBedpeOrdering(const BamAlignment &bam1,
                               BedpeMateFields &first,
                               BedpeMateFields &second,
                               bool mate1First)
{
    if (!mate1First) {
        if (ShouldSwapBedpeByCoordinate(first, second)) {
            swap(first, second);
        }
    }
    else {
        if (!bam1.IsFirstMate()) {
            swap(first, second);
        }
    }
}

static uint16_t CountMappedBedpeEnds(const BamAlignment &bam1,
                                     const BamAlignment &bam2)
{
    const unsigned mappedMask =
        static_cast<unsigned>(bam1.IsMapped()) |
        (static_cast<unsigned>(bam2.IsMapped()) << 1);
    unsigned count = 0;
    for (unsigned bit = 0; bit < 2; ++bit) {
        count += (mappedMask >> bit) & 1u;
    }
    return static_cast<uint16_t>(count);
}

static uint16_t BuildBedpeScore(const BamAlignment &bam1,
                                const BamAlignment &bam2,
                                const BedpeMateFields &first,
                                const BedpeMateFields &second,
                                bool useEditDistance)
{
    if (useEditDistance == false) {
        if (CountMappedBedpeEnds(bam1, bam2) == 2) {
            return min(bam1.MapQuality, bam2.MapQuality);
        }
        return 0;
    }

    uint32_t totalEditDistance = 0;
    if (bam1.IsMapped() == true && bam2.IsMapped() == true) {
        totalEditDistance = first.editDistance + second.editDistance;
    }
    else if (bam1.IsMapped() == true) {
        totalEditDistance = first.editDistance;
    }
    else if (bam2.IsMapped() == true) {
        totalEditDistance = second.editDistance;
    }
    return static_cast<uint16_t>(totalEditDistance);
}

static BedpeRecordContext BuildBedpeRecordContext(const BamAlignment &bam1,
                                                  const BamAlignment &bam2,
                                                  const RefVector &refs,
                                                  bool useEditDistance,
                                                  bool mate1First)
{
    BedpeRecordContext record;
    record.first = BuildBedpeMateFields(bam1, refs, useEditDistance);
    record.second = BuildBedpeMateFields(bam2, refs, useEditDistance);
    record.name = bam1.Name;
    record.score = BuildBedpeScore(bam1, bam2, record.first, record.second,
                                   useEditDistance);
    ApplyBedpeOrdering(bam1, record.first, record.second, mate1First);
    return record;
}

static vector<string> BuildBedpeFields(const BedpeRecordContext &record)
{
    vector<string> fields;
    fields.push_back(BuildBedValueText(record.first.chrom));
    fields.push_back(BuildBedValueText(record.first.start));
    fields.push_back(BuildBedValueText(record.first.end));
    fields.push_back(BuildBedValueText(record.second.chrom));
    fields.push_back(BuildBedValueText(record.second.start));
    fields.push_back(BuildBedValueText(record.second.end));
    fields.push_back(BuildBedValueText(record.name));
    fields.push_back(BuildBedValueText(record.score));
    fields.push_back(BuildBedValueText(record.first.strand));
    fields.push_back(BuildBedValueText(record.second.strand));
    return fields;
}

static string BuildBedpeRecord(const BedpeRecordContext &record)
{
    const vector<string> fields = BuildBedpeFields(record);
    BedFieldPolicy policy = BuildBedFieldPolicy(fields.size(), false, false);
    return BuildDelimitedBedRecord(fields, policy);
}

static vector<string> SplitBed12Color(const string &color)
{
    vector<string> channels;
    string channel;
    for (size_t i = 0; i < color.size(); ++i) {
        if (color[i] == ',') {
            channels.push_back(channel);
            channel.clear();
        }
        else {
            channel += color[i];
        }
    }
    channels.push_back(channel);
    return channels;
}

static bool ParseBed12ColorChannel(const string &text, int &value)
{
    if (text.empty()) {
        return false;
    }

    int parsed = 0;
    for (size_t i = 0; i < text.size(); ++i) {
        if (text[i] < '0' || text[i] > '9') {
            return false;
        }
        parsed = (parsed * 10) + (text[i] - '0');
    }
    if (parsed < 0 || parsed > 255) {
        return false;
    }
    value = parsed;
    return true;
}

static string NormalizeBed12Color(const string &color)
{
    vector<string> channels = SplitBed12Color(color);
    if (channels.size() != 3) {
        return color;
    }

    int red = 0;
    int green = 0;
    int blue = 0;
    if (!ParseBed12ColorChannel(channels[0], red) ||
        !ParseBed12ColorChannel(channels[1], green) ||
        !ParseBed12ColorChannel(channels[2], blue)) {
        return color;
    }

    // Preserve the caller's exact RGB text while keeping the validation path local
    // to BED12 formatting.
    return color;
}

static string BuildBed12Record(const Bed12RecordContext &record,
                               const vector<BED> &bedBlocks,
                               const BamAlignment &bam)
{
    string blockSizes = BuildBed12BlockField(bedBlocks, bam, record,
                                             BED12_BLOCK_FIELD_SIZE);
    string blockStarts = BuildBed12BlockField(bedBlocks, bam, record,
                                              BED12_BLOCK_FIELD_START);

    stringstream out;
    out << record.chrom << "\t"
        << record.chromStart << "\t"
        << record.chromEnd << "\t"
        << record.name << "\t"
        << record.score << "\t"
        << record.strand << "\t"
        << record.chromStart << "\t"
        << record.chromEnd << "\t"
        << record.color << "\t"
        << (int)bedBlocks.size() << "\t"
        << blockSizes << "\t"
        << blockStarts;
    return out.str();
}

string PrintTag(const BamAlignment &bam, const string &tag)
{
    uint32_t uTagValue;
    int32_t sTagValue;
    ostringstream value;
    if (bam.GetTag(tag, uTagValue))
        value << uTagValue;
    else if (bam.GetTag(tag, sTagValue))
        value << sTagValue;
    else {
        cerr << "The requested tag (" 
             << tag 
             << ") was not found in the BAM file.  Exiting\n";
        exit(1);
    }
    return value.str();
}

void PrintBed(const BamAlignment &bam,  const RefVector &refs, 
              bool useEditDistance, const string &bamTag, 
              bool obeySplits, bool splitOnDeletions,
              bool useCigar, bool useNovoalign, 
              bool useBWA) 
{
    BedRecordContext record = BuildBedRecordContext(bam, refs, bamTag,
                                                    useCigar, obeySplits);

    // report the entire BAM footprint as a single BED entry
    if (obeySplits == false) {
        
        if (!useNovoalign && !useBWA) {
            // report the alignment in BED6 format.
            cout << BuildStandardBedRecord(record);
        }
        else if (useNovoalign && !useBWA) {
            // special BED format for Hydra using Novoalign.
            cout << BuildHydraBedRecord(record, bam, BED_OUTPUT_NOVOALIGN);
        }
        else if (!useNovoalign && useBWA) {
            // special BED format for Hydra using BWA.
            cout << BuildHydraBedRecord(record, bam, BED_OUTPUT_BWA);
        }
    }
    // Report each chunk of the BAM alignment as a discrete BED entry
    // For example 10M100N10M would be reported as two seprate BED entries of
    // length 10
    else {
        // parse the CIGAR string and figure out the alignment blocks
        vector<BED> bedBlocks = BuildBamBedBlocks(bam, record.chrom,
                                                  splitOnDeletions);
        for (size_t i = 0; i < bedBlocks.size(); ++i) {
            cout << BuildSplitBedRecord(record, bedBlocks[i], bamTag, bam);
        }
    }
}


void PrintBed12(const BamAlignment &bam, const RefVector &refs, 
                bool useEditDistance, const string &bamTag, 
                bool obeySplits, bool splitOnDeletions,
                string color) 
{

    Bed12RecordContext record;
    record.chrom = refs.at(bam.RefID).RefName;
    record.name = BuildBamFeatureName(bam);
    record.score = BuildBed12Score(bam, bamTag);
    record.strand = bam.IsReverseStrand() ? "-" : "+";
    record.color = NormalizeBed12Color(color);
    record.chromStart = (CHRPOS)bam.Position;
    record.chromEnd = bam.GetEndPosition();
    record.splitOnDeletions = splitOnDeletions;
    record.cigarSummary = SummarizeBed12Cigar(bam.CigarData, splitOnDeletions);

    // parse the CIGAR string and figure out the alignment blocks
    vector<BED> bedBlocks = BuildBamBedBlocks(bam, record.chrom,
                                              splitOnDeletions);

    cout << BuildBed12Record(record, bedBlocks, bam) << endl;
}


void PrintBedPE(const BamAlignment &bam1, const BamAlignment &bam2, 
                const RefVector &refs, bool useEditDistance, bool mate1First) {

    BedpeRecordContext record = BuildBedpeRecordContext(bam1, bam2, refs,
                                                        useEditDistance,
                                                        mate1First);
    cout << BuildBedpeRecord(record);
}


// deprecated.
bool bamtobed_IsCorrectMappingForBEDPE (const BamAlignment &bam) {

    if ( (bam.RefID == bam.MateRefID) && (bam.InsertSize > 0) ) {
        return true;
    }
    else if ( (bam.RefID == bam.MateRefID) && 
              (bam.InsertSize == 0) && 
              bam.IsFirstMate() ) 
    {
        return true;
    }
    else if ( (bam.RefID != bam.MateRefID) && bam.IsFirstMate() ) {
        return true;
    }
    else return false;
}
