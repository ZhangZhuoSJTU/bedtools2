/*****************************************************************************
  BlockedIntervals.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#include "BlockedIntervals.h"

enum CigarBlockClass {
    CIGAR_BLOCK_ALIGNED,
    CIGAR_BLOCK_DELETION,
    CIGAR_BLOCK_SKIP,
    CIGAR_BLOCK_IGNORED,
    CIGAR_BLOCK_INVALID
};

static CigarBlockClass GetCigarBlockClass(char opType)
{
    switch (opType) {
        case 'M': case 'X': case '=':
            return CIGAR_BLOCK_ALIGNED;
        case 'D':
            return CIGAR_BLOCK_DELETION;
        case 'N':
            return CIGAR_BLOCK_SKIP;
        case 'I': case 'S': case 'P': case 'H':
            return CIGAR_BLOCK_IGNORED;
        default:
            return CIGAR_BLOCK_INVALID;
    }
}

static bool ShouldSplitCigarBlock(CigarBlockClass opClass,
                                  bool breakOnDeletionOps,
                                  bool breakOnSkipOps)
{
    return (opClass == CIGAR_BLOCK_DELETION && breakOnDeletionOps) ||
           (opClass == CIGAR_BLOCK_SKIP && breakOnSkipOps);
}

static void AppendBamBlock(const string &chrom,
                           const BamAlignment &bam,
                           const string &strand,
                           bedVector &bedBlocks,
                           CHRPOS start,
                           CHRPOS blockLength)
{
    bedBlocks.push_back( BED(chrom, start, start + blockLength,
                              bam.Name, ToString(bam.MapQuality), strand) );
}


void GetBamBlocks(const BamAlignment &bam,
                  const string &chrom,
                  bedVector &bedBlocks,
                  bool breakOnDeletionOps,
                  bool breakOnSkipOps)
{
    string strand;
    bam.IsReverseStrand() ? strand = "-" : strand = "+";
    CHRPOS currPosition = bam.Position;
    CHRPOS blockLength = 0;

    //  Rip through the CIGAR ops and figure out if there is more
    //  than one block for this alignment
    vector<CigarOp>::const_iterator cigItr = bam.CigarData.begin();
    vector<CigarOp>::const_iterator cigEnd = bam.CigarData.end();
    for (; cigItr != cigEnd; ++cigItr) {
        CigarBlockClass opClass = GetCigarBlockClass(cigItr->Type);
        if (opClass == CIGAR_BLOCK_ALIGNED ||
            ((opClass == CIGAR_BLOCK_DELETION || opClass == CIGAR_BLOCK_SKIP) &&
             !ShouldSplitCigarBlock(opClass, breakOnDeletionOps, breakOnSkipOps))) {
            blockLength += cigItr->Length;
        }
        else if (ShouldSplitCigarBlock(opClass, breakOnDeletionOps, breakOnSkipOps)) {
            if (blockLength > 0) {
                AppendBamBlock(chrom, bam, strand, bedBlocks, currPosition, blockLength);
                currPosition += blockLength;
                blockLength = 0;
            }
            currPosition += cigItr->Length;
        }
        else if (opClass == CIGAR_BLOCK_INVALID) {
            fprintf(stderr,"ERROR: Invalid Cigar op type \'%c\'.\n",cigItr->Type);   // shouldn't get here
            exit(1);
        }
    }
    AppendBamBlock(chrom, bam, strand, bedBlocks,
                   currPosition, blockLength);
}


void GetBedBlocks(const BED &bed, bedVector &bedBlocks) {

    // nothing to do if it is not a BED12 intervals
    if (bed.fields.size() != 12) {
        bedBlocks.push_back(bed);
        return;
    }

    int blockCount = atoi(bed.fields[9].c_str());
    if ( blockCount <= 0 ) {
        cerr << "Input error: found interval having <= 0 blocks." << endl;
        exit(1);
    }
    else {
        // get the comma-delimited strings for the BED12 block starts and block ends.
        string blockSizes(bed.fields[10]);
        string blockStarts(bed.fields[11]);

        vector<CHRPOS> sizes;
        vector<CHRPOS> starts;
        Tokenize(blockSizes, sizes, ',');
        Tokenize(blockStarts, starts, ',');

        if ( sizes.size() != (size_t) blockCount || starts.size() != (size_t) blockCount ) {
            cerr << "Input error: found interval with block-counts not matching starts/sizes on line." << endl;
            exit(1);
        }

        // add each BED block to the bedBlocks vector
        for (UINT i = 0; i < (UINT) blockCount; ++i) {
            CHRPOS blockStart = bed.start + starts[i];
            CHRPOS blockEnd   = bed.start + starts[i] + sizes[i];
            BED currBedBlock(bed.chrom, blockStart, blockEnd, 
                             bed.name, bed.score, bed.strand, bed.fields, bed.other_idxs);
            bedBlocks.push_back(currBedBlock);
        }
    }
}


int GetTotalBlockLength(const bedVector &bedBlocks) {
    int total_size = 0;
    bedVector::const_iterator blockItr = bedBlocks.begin();
    bedVector::const_iterator blockEnd = bedBlocks.end();
    for (; blockItr != blockEnd; ++blockItr) {
        total_size += blockItr->end - blockItr->start;
    }
    return total_size;
}
