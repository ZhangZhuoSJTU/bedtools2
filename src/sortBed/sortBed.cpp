/*****************************************************************************
  sortBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include <set>
#include "lineFileUtilities.h"
#include "sortBed.h"

using namespace SortBedConstants;

//
// Constructor: load the BED file into an (already start-sorted) per-chrom
// map and optionally echo the header before any sorted output.
//
BedSort::BedSort(string &bedFile, bool printHeader,string &faidxFile):_faidxFile(faidxFile) {
    _bedFile = bedFile;
    _bed = new BedFile(bedFile);
    
    _bed->loadBedFileIntoMapNoBin();
    // report the header first if asked.
    if (printHeader == true) {
        _bed->PrintHeader();
    }
}

//
// Destructor
//
BedSort::~BedSort(void) {
}


// Default sort: chromosome (ascending, via the map ordering) then start
// position (the per-chrom lists are already start-sorted on load).
void BedSort::SortBed() {

    // loop through each chromosome and merge their BED entries
    for (masterBedMapNoBin::iterator m = _bed->bedMapNoBin.begin(); m != _bed->bedMapNoBin.end(); ++m) {

        // bedList is already sorted by start position.
        const vector<BED> &bedList = m->second;

        for (unsigned int i = 0; i < bedList.size(); ++i) {
            _bed->reportBedNewLine(bedList[i]);
        }
    }
}

// Sort using an explicit chromosome order taken from a faidx/genome index
// file. Every chromosome present in the BED must be declared in the index,
// and the index must not declare any chromosome twice.
void BedSort::SortBedOnFaidx()
		{
		set<string> all_chromosomes;
		if(_faidxFile.empty())
				{
				cerr << "[sortBed] File for fasta index undefined." << endl;
				exit(EXIT_FAILURE);
				}
		/* read FAIDX file */
		ifstream faidx(_faidxFile.c_str(),ios::in);
		if(!faidx.is_open()) {
				cerr << "Cannot open \""<< _faidxFile << "\""
						<< strerror(errno)
						<< endl;
				exit(EXIT_FAILURE);
				}
		string line;
		while(getline(faidx,line,'\n')) {
				if(line.empty()) continue;
				string::size_type tab= line.find('\t');
				if(tab!=string::npos) {
						line.erase(tab);
						}
				if(all_chromosomes.find(line) != all_chromosomes.end()) {
						cerr  << "Chromosome \"" << line
									<<"\" defined twice in "
									<< _faidxFile
									<< endl;
						exit(EXIT_FAILURE);
						}
				_tid2chrom[_tid2chrom.size()]=line;
				all_chromosomes.insert(line);
				}
		faidx.close();
		/** end read FAIDX */
		
		//check BED chromosomes
		for (masterBedMapNoBin::iterator m = _bed->bedMapNoBin.begin(); m != _bed->bedMapNoBin.end(); ++m) {
					if( all_chromosomes.find(m->first) ==  all_chromosomes.end()) {
						cerr  << "Chromosome \"" << m->first
									<<"\" undefined in "
									<< _faidxFile
									<< endl;
						exit(EXIT_FAILURE);
						}
				}	
		
		
		//loop over each chromosome using the faidx index
		for(size_t tid=0; tid <_tid2chrom.size(); ++tid)
			{
			string chrom = _tid2chrom[tid];
			masterBedMapNoBin::iterator m = _bed->bedMapNoBin.find(chrom);
			
			if( m == _bed->bedMapNoBin.end() ) continue; //this chromosome is not present in BED
			
		    // bedList is already sorted by start position.
		    const vector<BED> &bedList = m->second;

		    for (unsigned int i = 0; i < bedList.size(); ++i) {
		        _bed->reportBedNewLine(bedList[i]);
			  }
			}
		
		}


// Sort every feature in the file by size (ascending), ignoring chromosome.
void BedSort::SortBedBySizeAsc() {

    vector<BED> masterList;
    masterList.reserve(MASTER_LIST_RESERVE);

    // loop through each chromosome and merge their BED entries
    for (masterBedMapNoBin::iterator m = _bed->bedMapNoBin.begin(); m != _bed->bedMapNoBin.end(); ++m) {

        // bedList is already sorted by start position.
        const vector<BED> &bedList = m->second;

        // add the entries from this chromosome to the current list
        for (unsigned int i = 0; i < bedList.size(); ++i) {
            masterList.push_back(bedList[i]);
        }
    }

    // sort the master list by size (asc.)
    sort(masterList.begin(), masterList.end(), sortBySizeAsc);

    // report the entries in ascending order
    for (unsigned int i = 0; i < masterList.size(); ++i) {
        _bed->reportBedNewLine(masterList[i]);
    }
}


// Sort every feature in the file by size (descending), ignoring chromosome.
void BedSort::SortBedBySizeDesc() {

    vector<BED> masterList;
    masterList.reserve(MASTER_LIST_RESERVE);

    // loop through each chromosome and merge their BED entries
    for (masterBedMapNoBin::iterator m = _bed->bedMapNoBin.begin(); m != _bed->bedMapNoBin.end(); ++m) {

        // bedList is already sorted by start position.
        const vector<BED> &bedList = m->second;

        // add the entries from this chromosome to the current list
        for (unsigned int i = 0; i < bedList.size(); ++i) {
            masterList.push_back(bedList[i]);
        }
    }

    // sort the master list by size (asc.)
    sort(masterList.begin(), masterList.end(), sortBySizeDesc);

    // report the entries in ascending order
    for (unsigned int i = 0; i < masterList.size(); ++i) {
        _bed->reportBedNewLine(masterList[i]);
    }
}

// Within each chromosome, sort the features by size (ascending).
void BedSort::SortBedByChromThenSizeAsc() {

    // loop through each chromosome and merge their BED entries
    for (masterBedMapNoBin::iterator m = _bed->bedMapNoBin.begin(); m != _bed->bedMapNoBin.end(); ++m) {

        // bedList is already sorted by start position.
        vector<BED> &bedList = m->second;
        sort(bedList.begin(), bedList.end(), sortBySizeAsc);

        for (unsigned int i = 0; i < bedList.size(); ++i) {
            _bed->reportBedNewLine(bedList[i]);
        }
    }
}


// Within each chromosome, sort the features by size (descending).
void BedSort::SortBedByChromThenSizeDesc() {

    // loop through each chromosome and merge their BED entries
    for (masterBedMapNoBin::iterator m = _bed->bedMapNoBin.begin(); m != _bed->bedMapNoBin.end(); ++m) {

        // bedList is already sorted by start position.
        vector<BED> &bedList = m->second;

        sort(bedList.begin(), bedList.end(), sortBySizeDesc);

        for (unsigned int i = 0; i < bedList.size(); ++i) {
            _bed->reportBedNewLine(bedList[i]);
        }
    }
}


// Within each chromosome, sort the features by score (ascending).
// Requires at least BED5 so that a score column exists.
void BedSort::SortBedByChromThenScoreAsc() {

    if (_bed->bedType >= MIN_BEDTYPE_FOR_SCORE) {
        // loop through each chromosome and merge their BED entries
        for (masterBedMapNoBin::iterator m = _bed->bedMapNoBin.begin(); m != _bed->bedMapNoBin.end(); ++m) {

            // bedList is already sorted by start position.
            vector<BED> &bedList = m->second;
            sort(bedList.begin(), bedList.end(), sortByScoreAsc);

            for (unsigned int i = 0; i < bedList.size(); ++i) {
                _bed->reportBedNewLine(bedList[i]);
            }
        }
    }
    else {
        cerr << "Error: Requested a sort by score, but your BED file does not appear to be in BED 5 format or greater.  Exiting." << endl;
        exit(EXIT_FAILURE_STATUS);
    }
}


// Within each chromosome, sort the features by score (descending).
// Requires at least BED5 so that a score column exists.
void BedSort::SortBedByChromThenScoreDesc() {

    if (_bed->bedType >= MIN_BEDTYPE_FOR_SCORE) {
        // loop through each chromosome and merge their BED entries
        for (masterBedMapNoBin::iterator m = _bed->bedMapNoBin.begin(); m != _bed->bedMapNoBin.end(); ++m) {

            // bedList is already sorted by start position.
            vector<BED> &bedList = m->second;
            sort(bedList.begin(), bedList.end(), sortByScoreDesc);

            for (unsigned int i = 0; i < bedList.size(); ++i) {
                _bed->reportBedNewLine(bedList[i]);
            }
        }
    }
    else {
        cerr << "Error: Requested a sort by score, but your BED file does not appear to be in BED 5 format or greater.  Exiting." << endl;
        exit(EXIT_FAILURE_STATUS);
    }
}

