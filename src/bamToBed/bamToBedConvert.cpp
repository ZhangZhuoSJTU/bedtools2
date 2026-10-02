/*****************************************************************************
  bamToBedConvert.cpp

  BAM to BED/BED12/BEDPE conversion routines, split out of bamToBed.cpp
******************************************************************************/
#include "bamToBed.h"

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
        switch (cigar[i].Type) {
            case ('M') :
            case ('=') :
            case ('X') :
            case ('I') :
            case ('D') :
            case ('N') :
            case ('S') :
            case ('H') :
            case ('P') :
                cigarString << cigar[i].Length << cigar[i].Type;
        }
    }
    return cigarString.str();
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
    // set the strand
    string strand = "+";
    if (bam.IsReverseStrand() == true) strand = "-";

    // set the name of the feature based on the sequence
    string name = bam.Name;
    if (bam.IsFirstMate() == true)  name += "/1";
    if (bam.IsSecondMate() == true) name += "/2";

    // get the unpadded (parm = false) end position based on the CIGAR
    unsigned int alignmentEnd = bam.GetEndPosition(false, false);

    // report the entire BAM footprint as a single BED entry
    if (obeySplits == false) {
        
        if (!useNovoalign && !useBWA) {
            // report the alignment in BED6 format.
            if (bamTag == "") {
                cout << refs.at(bam.RefID).RefName << "\t"
                     << bam.Position << "\t"
                     << alignmentEnd << "\t"
                     << name << "\t"
                     << bam.MapQuality << "\t"
                     << strand;
            }
            else {
                cout << refs.at(bam.RefID).RefName << "\t"
                     << bam.Position << "\t"
                     << alignmentEnd << "\t"
                     << name << "\t"
                     << PrintTag(bam, bamTag) << "\t"
                     << strand;
            }
            // does the user want CIGAR as well?
            if (useCigar == false) {
                printf("\n");
            }
            else {
                string cigar = BuildCigarString(bam.CigarData);
                printf("\t%s\n", cigar.c_str());
            }
        }
        else if (useNovoalign && !useBWA) {
            // special BED format for Hydra using Novoalign.
            uint32_t numMappings;

            if (!bam.GetTag("ZN", numMappings)) {
                // if ZN is missing, this means just one alignment was found.
                numMappings = 1;
            }
            cout << refs.at(bam.RefID).RefName << "\t"
                 << bam.Position << "\t"
                 << alignmentEnd << "\t"
                 << name << "\t"
                 << bam.MapQuality << "\t"
                 << PrintTag(bam, "NM") << "\t"
                 << numMappings << "\t"
                 << strand
                 << endl;
        }
        else if (!useNovoalign && useBWA) {
            // special BED format for Hydra using Novoalign.
            uint32_t numMappings;
            uint32_t x0 = 1;
            uint32_t x1 = 0;

            if (!bam.GetTag("X0", x0) && !bam.GetTag("X1", x1)) {
                // if ZN is missing, this means just one alignment was found.
                numMappings = 1;
            }
            else {
                numMappings = x0 + x1;
            }
            cout << refs.at(bam.RefID).RefName << "\t"
                 << bam.Position << "\t"
                 << alignmentEnd << "\t"
                 << name << "\t"
                 << bam.MapQuality << "\t"
                 << PrintTag(bam, "NM") << "\t"
                 << numMappings << "\t"
                 << strand
                 << endl;
        }
    }
    // Report each chunk of the BAM alignment as a discrete BED entry
    // For example 10M100N10M would be reported as two seprate BED entries of
    // length 10
    else {
        // parse the CIGAR string and figure out the alignment blocks
        vector<BED> bedBlocks;
        string chrom = refs.at(bam.RefID).RefName;
        // extract the block starts and lengths from the CIGAR string
        if (!splitOnDeletions)
            GetBamBlocks(bam, chrom, bedBlocks, false, true);
        else
            GetBamBlocks(bam, chrom, bedBlocks, true, true);


        unsigned int i;
        for (i = 0; i < bedBlocks.size(); ++i) {
            BED curr = bedBlocks[i];

            if (bamTag == "") {
                printf("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%d\t%s\n",
                       chrom.c_str(), 
                       curr.start,
                       curr.end,
                       name.c_str(),
                       (uint16_t)bam.MapQuality,
                       strand.c_str());
            }
            else {
                cout << chrom << "\t"
                     << bam.Position << "\t"
                     << curr.start << "\t"
                     << curr.end << "\t"
                     << name << "\t"
                     << PrintTag(bam, bamTag) << "\t"
                     << strand
                     << endl;            
            }
        }
    }
}


void PrintBed12(const BamAlignment &bam, const RefVector &refs, 
                bool useEditDistance, const string &bamTag, 
                bool obeySplits, bool splitOnDeletions,
                string color) 
{

    // set the strand
    string strand = "+";
    if (bam.IsReverseStrand()) strand = "-";

    // set the name of the feature based on the sequence
    string name = bam.Name;
    if (bam.IsFirstMate()) name += "/1";
    if (bam.IsSecondMate()) name += "/2";

    // parse the CIGAR string and figure out the alignment blocks
    vector<BED> bedBlocks;
    string chrom = refs.at(bam.RefID).RefName;
    CHRPOS alignmentEnd = bam.GetEndPosition();
    // extract the block starts and lengths from the CIGAR string
    if (!splitOnDeletions)
        GetBamBlocks(bam, chrom, bedBlocks, false, true);
    else
        GetBamBlocks(bam, chrom, bedBlocks, true, true);

    if (bamTag == "") {
        cout << refs.at(bam.RefID).RefName << "\t"
             << bam.Position << "\t"
             << alignmentEnd << "\t"
             << name << "\t"
             << bam.MapQuality << "\t"
             << strand << "\t";
    }
    else {
        cout << refs.at(bam.RefID).RefName << "\t"
             << bam.Position << "\t"
             << alignmentEnd << "\t"
             << name << "\t"
             << PrintTag(bam, bamTag) << "\t"
             << strand << "\t";
    }

    // write the colors, etc.
    printf("%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%d\t",
                               (CHRPOS)bam.Position, alignmentEnd,
                               color.c_str(), (int) bedBlocks.size());

    // now write the lengths portion
    unsigned int b;
    for (b = 0; b < bedBlocks.size() - 1; ++b) {
        printf("%" PRId_CHRPOS ",", bedBlocks[b].end - bedBlocks[b].start);
    }
    printf("%" PRId_CHRPOS "\t", bedBlocks[b].end - bedBlocks[b].start);

    // now write the starts portion
    for (b = 0; b < bedBlocks.size() - 1; ++b) {
        printf("%" PRId_CHRPOS ",", bedBlocks[b].start - bam.Position);
    }
    printf("%" PRId_CHRPOS "\n", bedBlocks[b].start - bam.Position);
}


void PrintBedPE(const BamAlignment &bam1, const BamAlignment &bam2, 
                const RefVector &refs, bool useEditDistance, bool mate1First) {

    // initialize BEDPE variables
    string chrom1, chrom2, strand1, strand2;
    int start1, start2, end1, end2;
    uint32_t editDistance1, editDistance2;
    start1 = start2 = end1 = end2 = -1;
    chrom1 = chrom2 = strand1 = strand2 = ".";
    editDistance1 = editDistance2 = 0;
    uint16_t minMapQuality = 0;

    // extract relevant info for end 1
    if (bam1.IsMapped()) {
        chrom1 = refs.at(bam1.RefID).RefName;
        start1 = bam1.Position;
        end1   = bam1.GetEndPosition(false);
        strand1 = "+";
        if (bam1.IsReverseStrand()) strand1 = "-";

        // extract the edit distance from the NM tag
        // if possible. otherwise, complain.
        if (useEditDistance == true && 
            bam1.GetTag("NM", editDistance1) == false) 
        {
            cerr << "The edit distance tag (NM) was not found in the BAM file.  Please disable -ed.  Exiting\n";
            exit(1);
        }
    }

    // extract relevant info for end 2
    if (bam2.IsMapped()) {
        chrom2 = refs.at(bam2.RefID).RefName;
        start2 = bam2.Position;
        end2   = bam2.GetEndPosition(false);
        strand2 = "+";
        if (bam2.IsReverseStrand()) strand2 = "-";

        // extract the edit distance from the NM tag
        // if possible. otherwise, complain.
        if (useEditDistance == true && 
            bam2.GetTag("NM", editDistance2) == false) 
        {
            cerr << "The edit distance tag (NM) was not found in the BAM file.  Please disable -ed.  Exiting\n";
            exit(1);
        }
    }

    // swap the ends if necessary
    if (!mate1First) 
    {
        if ( chrom1 > chrom2 || ((chrom1 == chrom2) && (start1 > start2)) ) {
            swap(chrom1, chrom2);
            swap(start1, start2);
            swap(end1, end2);
            swap(strand1, strand2);
        }
    }
    // if -mate1First, make sure the first block is always mate 1
    else 
    {
        if (!bam1.IsFirstMate()) {
            swap(chrom1, chrom2);
            swap(start1, start2);
            swap(end1, end2);
            swap(strand1, strand2);
        }
    }

    // report BEDPE using min mapQuality
    if (useEditDistance == false) {
        // compute the minimum mapping quality b/w the two ends of the pair.
        if (bam1.IsMapped() == true && bam2.IsMapped() == true)
            minMapQuality = min(bam1.MapQuality, bam2.MapQuality);

        printf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\n",
                chrom1.c_str(), start1, 
                end1, chrom2.c_str(), 
                start2, end2,
                bam1.Name.c_str(), minMapQuality, 
                strand1.c_str(), strand2.c_str());
    }
    // report BEDPE using total edit distance
    else {
        uint16_t totalEditDistance = 0;
        if (bam1.IsMapped() == true && bam2.IsMapped() == true)
            totalEditDistance = editDistance1 + editDistance2;
        else if (bam1.IsMapped() == true)
            totalEditDistance = editDistance1;
        else if (bam2.IsMapped() == true)
            totalEditDistance = editDistance2;

        printf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\n",
                chrom1.c_str(), start1, 
                end1, chrom2.c_str(), 
                start2, end2,
                bam1.Name.c_str(), totalEditDistance, 
                strand1.c_str(), strand2.c_str());
    }
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
