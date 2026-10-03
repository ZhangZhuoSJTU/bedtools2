/*****************************************************************************
  GenomeFile.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "GenomeFile.h"

// Process-exit status used when the genome file is missing or malformed.
static const int GENOMEFILE_ERROR_EXIT_CODE = 1;

// A genome-file line must have at least a chromosome name and a length.
static const size_t MIN_GENOME_FILE_FIELDS = 2;

// Column indices within a genome-file line.
static const size_t CHROM_NAME_FIELD = 0;
static const size_t CHROM_SIZE_FIELD = 1;

// Base used when parsing a chromosome length with strtol.
static const int CHROM_SIZE_NUMERIC_BASE = 10;

// Leading token that marks a line as a comment, to be skipped.
static const char * const COMMENT_PREFIX = "#";


GenomeFile::GenomeFile(const string &genomeFile) {
    _genomeFile = genomeFile;
    loadGenomeFileIntoMap();
}

GenomeFile::GenomeFile(const RefVector &genome) {
    for (size_t i = 0; i < genome.size(); ++i) {
        string chrom = genome[i].RefName;
        int length = genome[i].RefLength;
        
        _chromSizes[chrom] = length;
        _chromList.push_back(chrom);
    }
}

// Destructor
GenomeFile::~GenomeFile(void) {
}


// Read the genome file line by line, storing each valid chromosome's size and
// building the cumulative start-offset table used by projectOnGenome().
void GenomeFile::loadGenomeFileIntoMap() {

    string genomeLine;
    int lineNum = 0;
    vector<string> genomeFields;            // vector for a GENOME entry

    // open the GENOME file for reading
    ifstream genome(_genomeFile.c_str(), ios::in);
    if ( !genome ) {
        cerr << "Error: The requested genome file (" << _genomeFile << ") could not be opened. Exiting!" << endl;
        exit (GENOMEFILE_ERROR_EXIT_CODE);
    }

    while (getline(genome, genomeLine)) {

        Tokenize(genomeLine,genomeFields);  // load the fields into the vector
        lineNum++;

        // ignore a blank line
        if (genomeFields.size() > 0) {
            // Skip comment lines (those whose first field begins with '#').
            if (genomeFields[CHROM_NAME_FIELD].find(COMMENT_PREFIX) != 0) {
                // we need at least 2 columns
                if (genomeFields.size() >= MIN_GENOME_FILE_FIELDS) {
                    char *sizeParseEnd;
                    long parsedSize;
                    // make sure the second column is numeric.
                    parsedSize = strtol(genomeFields[CHROM_SIZE_FIELD].c_str(), &sizeParseEnd, CHROM_SIZE_NUMERIC_BASE);
                    (void)parsedSize; // only used to validate numeric-ness via sizeParseEnd

                    // strtol  will set sizeParseEnd to the start of the string if non-integral, base 10
                    if (sizeParseEnd != genomeFields[CHROM_SIZE_FIELD].c_str()) {
                        string chrom       = genomeFields[CHROM_NAME_FIELD];
                        CHRPOS size           = atol(genomeFields[CHROM_SIZE_FIELD].c_str());
                        _chromSizes[chrom] = size;
                        _chromList.push_back(chrom);
                        _startOffsets.push_back(_genomeLength);
                        _genomeLength += size;
                    }
                }
                else {
                    cerr << "Less than the req'd two fields were encountered in the genome file (" << _genomeFile << ")";
                    cerr << " at line " << lineNum << ".  Exiting." << endl;
                    exit (GENOMEFILE_ERROR_EXIT_CODE);
                }
            }
        }
        genomeFields.clear();
    }
}

// Translate a genome-wide linear position into the (chromosome name, offset
// within that chromosome) pair it corresponds to.
pair<string, CHRPOS> GenomeFile::projectOnGenome(CHRPOS genome_pos) {
    // search for the chrom that the position belongs on.
    // add 1 to genome position b/c of zero-based, half open.
    vector<CHRPOS>::const_iterator low =
        lower_bound(_startOffsets.begin(), _startOffsets.end(), genome_pos + 1);
    
    // use the iterator to identify the appropriate index 
    // into the chrom name and start vectors
    CHRPOS i = CHRPOS(low-_startOffsets.begin());
    string chrom = _chromList[i - 1];
    CHRPOS start = genome_pos - _startOffsets[i - 1];
    return make_pair(chrom, start);
}
    
// Return the size of the named chromosome, or (uint64_t)-1 if it is unknown.
uint64_t GenomeFile::getChromSize(const string &chrom) {
    chromToSizes::const_iterator chromIt = _chromSizes.find(chrom);
    if (chromIt != _chromSizes.end())
        return (uint64_t)_chromSizes[chrom];
    else
        return -1;  // chrom not found.
}

uint64_t GenomeFile::getGenomeSize(void) {
    return (uint64_t)_genomeLength;
}

vector<string> GenomeFile::getChromList() {
    return _chromList;
}

int GenomeFile::getNumberOfChroms() {
    return _chromList.size();
}

string GenomeFile::getGenomeFileName() {
    return _genomeFile;
}
