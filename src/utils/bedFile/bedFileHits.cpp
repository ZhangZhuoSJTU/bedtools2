/*****************************************************************************
  bedFileHits.cpp

  BedFile hit counting routines, split out of bedFile.cpp
******************************************************************************/
#include "bedFile.h"

void BedFile::countHits(const BED &a, bool sameStrand, bool diffStrand, bool countsOnly) {

    BIN startBin, endBin;
    startBin = (a.start >> _binFirstShift);
    endBin = ((a.end-1) >> _binFirstShift);

    // loop through each bin "level" in the binning hierarchy
    for (BINLEVEL i = 0; i < _binLevels; ++i) {

        // loop through each bin at this level of the hierarchy
        BIN offset = _binOffsetsExtended[i];
        for (BIN j = (startBin+offset); j <= (endBin+offset); ++j) {
            // loop through each feature in this chrom/bin and 
            // see if it overlaps with the feature that was passed in.  
            // if so, add the feature to the list of hits.
            vector<BEDCOV>::iterator bedItr = bedCovMap[a.chrom][j].begin();
            vector<BEDCOV>::iterator bedEnd = bedCovMap[a.chrom][j].end();
            for (; bedItr != bedEnd; ++bedItr) {
                
                bool strands_are_same = (a.strand == bedItr->strand);
                // skip the hit if not on the same strand (and we care)
                if ((sameStrand == true && strands_are_same == false) ||
                    (diffStrand == true && strands_are_same == true)
                   ) 
                {
                    continue;
                }
                else if (overlaps(bedItr->start, bedItr->end, a.start, a.end) 
                         > 0) 
                {
                    bedItr->count++;
                    if (countsOnly == false) {
                        if (a.zeroLength == false) {
                            bedItr->depthMap[a.start+1].starts++;
                            bedItr->depthMap[a.end].ends++;
                        }
                        else {
                            // correct for the fact that we artificially 
                            // expanded the zeroLength feature
                            bedItr->depthMap[a.start+2].starts++;
                            bedItr->depthMap[a.end-1].ends++;                   
                        }

                        if (a.start < bedItr->minOverlapStart) {
                            bedItr->minOverlapStart = a.start;
                        }
                    }
                }
            }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
}


void BedFile::countSplitHits(const vector<BED> &bedBlocks, bool sameStrand, bool diffStrand, bool countsOnly) {

    // set to track the distinct B features that had coverage.
    // we'll update the counts of coverage for these features by one
    // at the end of this function to avoid over-counting.
    set< vector<BEDCOV>::iterator > validHits;

    vector<BED>::const_iterator blockItr  = bedBlocks.begin();
    vector<BED>::const_iterator blockEnd  = bedBlocks.end();
    for (; blockItr != blockEnd; ++blockItr) {

        BIN startBin, endBin;
        startBin = (blockItr->start >> _binFirstShift);
        endBin = ((blockItr->end-1) >> _binFirstShift);

        // loop through each bin "level" in the binning hierarchy
        for (BINLEVEL i = 0; i < _binLevels; ++i) {

            // loop through each bin at this level of the hierarchy
            BIN offset = _binOffsetsExtended[i];
            for (BIN j = (startBin+offset); j <= (endBin+offset); ++j) {
                // loop through each feature in this chrom/bin and see if it 
                // overlaps with the feature that was passed in.  
                // if so, add the feature to the list of hits.
                vector<BEDCOV>::iterator 
                    bedItr = bedCovMap[blockItr->chrom][j].begin();
                vector<BEDCOV>::iterator 
                    bedEnd = bedCovMap[blockItr->chrom][j].end();
                for (; bedItr != bedEnd; ++bedItr) {
                    bool strands_are_same = 
                        (blockItr->strand == bedItr->strand);
                    // skip the hit if not on the same strand (and we care)
                    if ((sameStrand == true && strands_are_same == false) ||
                        (diffStrand == true && strands_are_same == true)
                       ) 
                    {
                        continue;
                    }
                    else if (overlaps(bedItr->start, bedItr->end, 
                                      blockItr->start, blockItr->end) > 0) 
                    {
                        if (countsOnly == false) {
                            if (blockItr->zeroLength == false) {
                                bedItr->depthMap[blockItr->start+1].starts++;
                                bedItr->depthMap[blockItr->end].ends++;
                            }
                            else {
                                // correct for the fact that we artificially 
                                // expanded the zeroLength feature
                                bedItr->depthMap[blockItr->start+2].starts++;
                                bedItr->depthMap[blockItr->end-1].ends++;
                            }
                        }

                        validHits.insert(bedItr);
                        if (blockItr->start < bedItr->minOverlapStart)
                            bedItr->minOverlapStart = blockItr->start;
                    }
                }
            }
            startBin >>= _binNextShift;
            endBin >>= _binNextShift;
        }
    }
    // incrment the count of overlap by one for each B feature that overlapped
    // the current passed hit.  This is necessary to prevent over-counting for
    // each "split"" of a single read.
    set< vector<BEDCOV>::iterator >::iterator validHitsItr = validHits.begin();
    set< vector<BEDCOV>::iterator >::iterator validHitsEnd = validHits.end();
    for (; validHitsItr != validHitsEnd; ++validHitsItr)
        // the validHitsItr points to another itr, hence 
        // the (*itr)-> dereferencing. ugly, but that's C++.
        (*validHitsItr)->count++;
}


void BedFile::countListHits(const BED &a, int index, bool sameStrand, bool diffStrand) {

    BIN startBin, endBin;
    startBin = (a.start >> _binFirstShift);
    endBin = ((a.end-1) >> _binFirstShift);

    // loop through each bin "level" in the binning hierarchy
    for (BINLEVEL i = 0; i < _binLevels; ++i) {

        // loop through each bin at this level of the hierarchy
        BIN offset = _binOffsetsExtended[i];
        for (BIN j = (startBin+offset); j <= (endBin+offset); ++j) {

            // loop through each feature in this chrom/bin and see if it 
            // overlaps with the feature that was passed in.  if so, 
            // add the feature tothe list of hits.
            vector<BEDCOVLIST>::iterator 
                bedItr = bedCovListMap[a.chrom][j].begin();
            vector<BEDCOVLIST>::iterator 
                bedEnd = bedCovListMap[a.chrom][j].end();
            for (; bedItr != bedEnd; ++bedItr) {

                bool strands_are_same = (a.strand == bedItr->strand);
                // skip the hit if not on the same strand (and we care)
                if ((sameStrand == true && strands_are_same == false) ||
                    (diffStrand == true && strands_are_same == true)
                   ) 
                {
                    continue;
                }
                else if (overlaps(bedItr->start, bedItr->end, 
                                  a.start, a.end) > 0) 
                {
                    bedItr->counts[index]++;
                    if (a.zeroLength == false) {
                        bedItr->depthMapList[index][a.start+1].starts++;
                        bedItr->depthMapList[index][a.end].ends++;
                    }
                    else {
                        // correct for the fact that we artificially expanded 
                        // the zeroLength feature
                        bedItr->depthMapList[index][a.start+2].starts++;
                        bedItr->depthMapList[index][a.end-1].ends++;
                    }

                    if (a.start < bedItr->minOverlapStarts[index]) {
                        bedItr->minOverlapStarts[index] = a.start;
                    }
                }
            }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
}
