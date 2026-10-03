
/*****************************************************************************
  NewGenomeFile.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef NEW_GENOMEFILE_H
#define NEW_GENOMEFILE_H

#include <algorithm> // for bsearch lower_bound()

#include "BedtoolsTypes.h"

#include "api/BamReader.h"
#include "api/BamAux.h"


// NewGenomeFile models a bedtools "genome" (a set of chromosomes and their
// lengths), loaded either from a two-column text file / .fai index or from the
// reference list of an open BAM file. Besides simple name->size lookups it
// maintains a cumulative offset table so that per-chromosome coordinates can be
// projected onto, and recovered from, a single linear genome-wide coordinate.
class NewGenomeFile {

public:

     NewGenomeFile(const string &genomeFileName);
     NewGenomeFile(const BamTools::RefVector &genome);
    ~NewGenomeFile(void);

    // load a GENOME file into a map keyed by chrom. value is a pair<int, int> of id and size.
    void loadGenomeFileIntoMap();
    
    bool projectOnGenome(CHRPOS genome_pos, string &chrom, CHRPOS &start);
    
    CHRPOS getGenomeSize(void) const { return _genomeLength; }                // return the total size of the genome
    CHRPOS getChromSize(const string &chrom);  // return the size of a chromosome
    CHRPOS getChromSize(const string &chrom) const;  // return the size of a chromosome
    uint32_t getChromId(const string &chrom); // return chromosome's sort order
    const vector<string> &getChromList() const { return _chromList; }  // return a list of chrom names
    uint32_t getNumberOfChroms() const { return (uint32_t)(_chromList.size() -1); }//the -1 excludes the blank chrom added for unmapped reads
    const string &getGenomeFileName() const { return _genomeFileName; }
    bool hasChrom(const string &chrom) const { return _chromSizeIds.find(chrom) != _chromSizeIds.end(); }




private:
    string  _genomeFileName;   // path to the genome file (or "-"/"stdin")
    istream   *_genomeFile;     // stream the genome is read from
    // Maps a chromosome name to a (size, sort-order id) pair.
    typedef map<string, pair<CHRPOS, int> > lookupType;
    lookupType _chromSizeIds;
    vector<string> _chromList;  // chromosome names in file order
    int _maxId;                 // highest sort-order id assigned so far

    // projecting chroms onto a single coordinate system
    CHRPOS _genomeLength;        // total length of all chromosomes combined
    vector<CHRPOS> _startOffsets;// cumulative start offset of each chromosome

    //cache members for quick lookup of the most recently queried chromosome
    string _currChromName;
    CHRPOS _currChromSize;
    int _currChromId;

};

#endif /* GENOMEFILE_H */
