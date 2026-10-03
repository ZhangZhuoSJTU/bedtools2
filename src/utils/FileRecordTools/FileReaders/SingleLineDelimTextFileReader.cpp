#include "SingleLineDelimTextFileReader.h"
#include <iostream>
#include <limits.h>
#include "BufferedStreamMgr.h"
#include "ParseTools.h"

bool abs_cmp(int i, int j) { return abs(i)<abs(j); }

SingleLineDelimTextFileReader::SingleLineDelimTextFileReader(int numFields, char delimChar)
: _numFields(numFields),
  _delimChar(delimChar),
  _fullHeaderFound(false),
  _currDataPos(0),
  _lineNum(0),
  _inheader(false)
{
	_delimPositions = new int[numFields +1];
}

SingleLineDelimTextFileReader::~SingleLineDelimTextFileReader()
{
	delete [] _delimPositions;
	_delimPositions = NULL;
}

bool SingleLineDelimTextFileReader::readEntry()
{
	if (!_isFileOpen) {
		return false;
	}

	if (_bufStreamMgr->eof()) {
		return false;
	}
	if (!_bufStreamMgr->getLine(_sLine)) {
		return false;
	}
	_lineNum++;
	if (_sLine.empty()) {
		return false;
	}

	//scan the whole header in one call.
	bool wasHeader = false;
	while (detectAndHandleHeader()) { //header line
		_lineNum++;
		if (!_bufStreamMgr->getLine(_sLine)) {
			return false;
		}
	}
	//after the first time we find a header, any other header line
	//doesn't count, and should not be added to the header, because it is a
	// "commented out record."
	_fullHeaderFound = true;

	// A line made up entirely of whitespace carries no record; skip it
	// rather than feeding it to the delimiter scanner.
	bool hasContent = false;
	int lineLength = (int)_sLine.size();
	for (int i = 0; i < lineLength; i++) {
		if (!isspace(_sLine[i])) {
			hasContent = true;
			break;
		}
	}
	if (!hasContent) {
		return false;
	}

	if (wasHeader) {
		return true;
	}

	// Build the field position table for this line; a field-count mismatch
	// is reported from within and terminates the run.
	if (!findDelimiters()) {
		return false;
	}
	return true;
}


int SingleLineDelimTextFileReader::startOffset(int fieldNum) const {
	// The table entry preceding a field holds its leading delimiter's
	// position; the field's own first character is one byte later.
	LineOffset delimPos = _delimPositions[fieldNum];
	return (int)delimPos + 1;
}

int SingleLineDelimTextFileReader::endOffset(int fieldNum) const {
	// The field ends where the next delimiter begins.
	LineOffset delimPos = _delimPositions[fieldNum + 1];
	return (int)delimPos;
}

int SingleLineDelimTextFileReader::fieldLength(int fieldNum) const {
	return endOffset(fieldNum) - startOffset(fieldNum);
}

void SingleLineDelimTextFileReader::getField(int fieldNum, string &str) const {
	str.assign(_sLine.c_str() + startOffset(fieldNum), fieldLength(fieldNum));
}


void SingleLineDelimTextFileReader::getField(int fieldNum, CHRPOS &val) const {
	string temp;
	getField(fieldNum, temp);
	val = str2chrPos(temp.c_str());
}

void SingleLineDelimTextFileReader::getField(int fieldNum, char &val) const {
	val = _sLine[startOffset(fieldNum)];
}

void SingleLineDelimTextFileReader::appendField(int fieldNum, string &str) const {
	str.append(_sLine.c_str() + startOffset(fieldNum), fieldLength(fieldNum));
}

bool SingleLineDelimTextFileReader::detectAndHandleHeader()
{
	//not sure why the linker is giving me a hard time about
	//passing a non-const string to isHeaderLine, but
	//this const ref is a workaround.
	const string &sLine2 = _sLine;
	if (!isHeaderLine(sLine2) && (!(_inheader && _lineNum==1))) {
		return false;
	}
	if (!_fullHeaderFound) {
		_header += _sLine;
		_header += "\n"; //add new line, since it was chomped by getline
	}
	return true;
}

bool SingleLineDelimTextFileReader::findDelimiters() {
	// Reset the position table for the new line before rescanning.
	memset(_delimPositions, 0, (_numFields + 1) * sizeof(int));

	// Walk the line once, recording the position of each field boundary.
	// Each table slot stores the index just before a field starts; the
	// trailing slot is closed off with the line length below.
	int fieldCount = 0;
	const char *lineStart = _sLine.c_str();
	const char *cursor = lineStart;
	const int delimByte = _delimChar;

	while (cursor) {
		_delimPositions[fieldCount++] = cursor - lineStart - 1;
		cursor = strchr(cursor, delimByte);
		if (cursor) {
			cursor++;
		}
	}
	_delimPositions[fieldCount] = _sLine.size();

	// A data line is only valid if it carries exactly the number of fields
	// the file's type was established with. Anything else (e.g. a stray
	// trailing tab) is a malformed record.
	if (fieldCount != _numFields) {
		cerr << "Error: line number "
			 << _lineNum << " of file "
			 << _filename << " has "
			 << fieldCount << " fields, but "
			 << _numFields
			 << " were expected."
			 << endl;
		exit(1);
	}
	return true;
}

CHRPOS SingleLineDelimTextFileReader::getVcfSVlen() {
	// tokenize the INFO field
	string info_str;
	vector<string> infofields;
	getField(VCF_TAG_FIELD, info_str);
	Tokenize(info_str, infofields, ';');

	// FOO=BAR;BIZ=BAM;SVLEN=100;END=200
	for (vector<string>::iterator f = infofields.begin(); f != infofields.end(); ++f) {
    	if (*f == ".") {
    		continue;
        }
        vector<string> keytoval;
        Tokenize(*f, keytoval, '=');
        //hey->val
        //SVLEN->100
        if (keytoval.size() == 2) {
        	if (keytoval.at(0) == "SVLEN") {
        		vector<CHRPOS> svlens;
        		Tokenize(keytoval.at(1), svlens, ',');
        		// are the multiple SVLENS?
        		if (svlens.size() == 1) {
        			return abs(svlens[0]);
        		}
        		else {
        			// return the abs_max SVLEN
        			CHRPOS max_len = *max_element(svlens.begin(),svlens.end(), abs_cmp);
        			return abs(max_len);
        		}
        	}
        	else if (keytoval.at(0) == "END") {
        		string start_str;
        		getField(1, start_str);
        		// length is END - POS + 1
        		return (int)(str2chrPos(keytoval.at(1)) - str2chrPos(start_str) + 1);
        	}
		}
    }
    // not found
    return INT_MIN;
}
