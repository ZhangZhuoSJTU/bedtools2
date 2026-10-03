/*****************************************************************************
  NewGenomeFile.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#include "NewGenomeFile.h"
#include "ParseTools.h"
#include "Tokenizer.h"
#include <fstream>

// Process-exit status used when the genome file cannot be opened or is
// malformed. Replaces the bare literal "1" used throughout this file.
static const int GENOMEFILE_ERROR_EXIT_CODE = 1;

// A bedtools genome file must supply at least a chromosome name and its
// length, i.e. two whitespace-delimited columns per line.
static const int MIN_GENOME_FILE_FIELDS = 2;

// Column indices within a genome-file line.
static const int CHROM_NAME_FIELD = 0;
static const int CHROM_SIZE_FIELD = 1;

// Chromosome lengths must be strictly positive; a length of zero is an error.
static const CHRPOS INVALID_CHROM_SIZE = 0;

// Initial value of _maxId before any chromosome has been seen. Each chromosome
// bumps _maxId, so a value of -1 after loading means the file was empty.
static const int NO_CHROMS_LOADED = -1;

// Name and size of the synthetic "unmapped" chromosome appended to every
// genome so that BAM records with no mapping position do not error out.
static const char * const UNMAPPED_CHROM_NAME = "";
static const CHRPOS UNMAPPED_CHROM_SIZE = 0;

NewGenomeFile::NewGenomeFile(const string &genomeFilename)
: _maxId(NO_CHROMS_LOADED)
{
    _genomeFileName = genomeFilename;
    loadGenomeFileIntoMap();
}

NewGenomeFile::NewGenomeFile(const BamTools::RefVector &refVector)
: _maxId(NO_CHROMS_LOADED)
{
	// Seed the chromosome lookup directly from the BAM reference list.
	size_t refIdx = 0;
    for (; refIdx < refVector.size(); ++refIdx) {
        string chrom = refVector[refIdx].RefName;
        CHRPOS length = refVector[refIdx].RefLength;
        _maxId++;
        _chromSizeIds[chrom] = pair<CHRPOS, CHRPOS>(length, _maxId);
		_chromList.push_back(chrom);
    }
	// Special: BAM files can have unmapped reads, which show as no chromosome, or an empty chrom string.
	// Add in an empty chrom so these don't error.
	_maxId++;
	_chromSizeIds[UNMAPPED_CHROM_NAME] = pair<CHRPOS, int>(UNMAPPED_CHROM_SIZE, _maxId);
	_chromList.push_back(UNMAPPED_CHROM_NAME);

}

// Destructor
NewGenomeFile::~NewGenomeFile(void) {
}


// Parse the genome file named by _genomeFileName and populate the chromosome
// lookup table, the ordered chromosome list and the cumulative start-offset
// table used to project chromosome coordinates onto a single linear genome.
void NewGenomeFile::loadGenomeFileIntoMap() {

	// Read from stdin when requested, otherwise open the named file.
	if (_genomeFileName == "-" || _genomeFileName == "stdin")
	{
		_genomeFile = &cin;
	}
	else
	{
		_genomeFile = new ifstream(_genomeFileName.c_str(), ios::in);
		if (!_genomeFile->good())
		{
			cerr << "Error: Can't open genome file " << _genomeFileName << ". Exiting..." << endl;
			exit(GENOMEFILE_ERROR_EXIT_CODE);
		}
	}

	string line;
	Tokenizer fieldTokens;
	CHRPOS chrSize = 0;
	_genomeLength = 0;
	string chrName;
	while (getline(*_genomeFile, line)) {
		trimNewlines(line);
		chrSize = 0;
		chrName.clear();
		int numFields = fieldTokens.tokenize(line.c_str());

		// allow use of .fai files: skip lines that lack the required columns.
		if (numFields < MIN_GENOME_FILE_FIELDS) {
			continue;
		}
		chrName = fieldTokens.getElem(CHROM_NAME_FIELD);
		chrSize = str2chrPos(fieldTokens.getElem(CHROM_SIZE_FIELD));
		if (chrSize == INVALID_CHROM_SIZE) {
			cerr << "Error: "
			     << chrName
			     << " has length equal to 0. Each chromosome must have non-zero length. Exiting."
			     << endl;
			exit(GENOMEFILE_ERROR_EXIT_CODE);
		}
		_maxId++;
		_chromSizeIds[chrName] = pair<CHRPOS, int>(chrSize, _maxId);
		// Record where this chromosome begins in the linearised genome, then
		// advance the running total by this chromosome's length.
		_startOffsets.push_back(_genomeLength);
		_genomeLength += chrSize;
		_chromList.push_back(chrName);
	}
	if (_maxId == NO_CHROMS_LOADED) {
		cerr << "Error: The genome file " << _genomeFileName << " has no valid entries (are you sure it's a 2-column bedtools genome file). Exiting." << endl;
		exit(GENOMEFILE_ERROR_EXIT_CODE);
	}
	// Special: BAM files can have unmapped reads, which show as no chromosome, or an empty chrom string.
	// Add in an empty chrom so these don't error.
	_maxId++;
	_chromSizeIds[UNMAPPED_CHROM_NAME] = pair<CHRPOS, int>(UNMAPPED_CHROM_SIZE, _maxId);
	_chromList.push_back(UNMAPPED_CHROM_NAME);


	_startOffsets.push_back(_genomeLength); //insert the final length as the last element
	//to help with the lower_bound call in the projectOnGenome method.
}

// Map a position in the linearised genome back to a (chromosome, offset)
// pair. Returns false (without writing the out-params) when genome_pos lies
// beyond the end of the genome.
bool NewGenomeFile::projectOnGenome(CHRPOS genome_pos, string &chrom, CHRPOS &start) {
    // search for the chrom that the position belongs on.
    // add 1 to genome position b/c of zero-based, half open.
    vector<CHRPOS>::const_iterator low =
        lower_bound(_startOffsets.begin(), _startOffsets.end(), genome_pos + 1);
    
    // use the iterator to identify the appropriate index 
    // into the chrom name and start vectors
    CHRPOS i = CHRPOS(low-_startOffsets.begin());
    if (i >= (CHRPOS)_chromList.size()) {
    	return false; //position not on genome
    }
    chrom = _chromList[i - 1];
    start = genome_pos - _startOffsets[i - 1];
    return true;
}
    
// Return the size of the named chromosome. This (non-const) overload also
// refreshes the single-entry lookup cache (_currChrom*) so that repeated
// queries for the same chromosome avoid the map lookup.
CHRPOS NewGenomeFile::getChromSize(const string &chrom) {
	if (chrom == _currChromName) {
		return _currChromSize;
	}
	lookupType::const_iterator iter= _chromSizeIds.find(chrom);
    if (iter != _chromSizeIds.end()) {
    	_currChromName = iter->first;
    	_currChromSize = iter->second.first;
    	_currChromId = iter->second.second;
    	return _currChromSize;
    }
    cerr << "Error: chromosome " << chrom << " is not in the genome file " << _genomeFileName << ". Exiting." << endl;
    return INT_MAX;
}

// const overload of getChromSize: performs the same lookup but cannot update
// the cache, so it only consults it.
CHRPOS NewGenomeFile::getChromSize(const string &chrom) const {
	if (chrom == _currChromName) {
		return _currChromSize;
	}
	lookupType::const_iterator iter= _chromSizeIds.find(chrom);
    if (iter != _chromSizeIds.end()) {
    	return iter->second.first;
    }
    cerr << "Error: chromosome " << chrom << " is not in the genome file " << _genomeFileName << ". Exiting." << endl;
    return INT_MAX;
}

// Return the sort-order id assigned to the named chromosome (the order in
// which it appeared in the genome file). Updates the lookup cache, and unlike
// getChromSize, aborts the program if the chromosome is unknown.
uint32_t NewGenomeFile::getChromId(const string &chrom) {
	if (chrom == _currChromName) {
		return _currChromId;
	}
	lookupType::const_iterator iter= _chromSizeIds.find(chrom);
    if (iter != _chromSizeIds.end()) {
    	_currChromName = iter->first;
    	_currChromSize = iter->second.first;
    	_currChromId = iter->second.second;
    	return _currChromId;
    } else {
    	cerr << "Error: requested chromosome " << chrom << " does not exist in the genome file " << _genomeFileName << ". Exiting." << endl;
    	exit(GENOMEFILE_ERROR_EXIT_CODE);
    }
    return INT_MAX;
}
