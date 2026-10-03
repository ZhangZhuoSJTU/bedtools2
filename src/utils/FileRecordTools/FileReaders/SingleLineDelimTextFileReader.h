/*
 * SingleLineTextFileReader.h
 *
 *  Created on: Nov 8, 2012
 *      Author: nek3d
 */

#ifndef SINGLELINETEXTFILEREADER_H_
#define SINGLELINETEXTFILEREADER_H_

#include <algorithm>
#include "FileReader.h"
#include "string.h"
#include "lineFileUtilities.h"

typedef int64_t CHRPOS;

// Byte offset of a field boundary within the line currently held by the
// reader. Records are scanned one line at a time, so an offset only ever
// addresses within that single line's buffer.
typedef short LineOffset;

class SingleLineDelimTextFileReader : public FileReader {
public:
	//Allow VCF records to access a specialized private method.
	//See end of class declaration for details.
	friend class VcfRecord;

	SingleLineDelimTextFileReader(int numFields, char delimChar = '\t');
	~SingleLineDelimTextFileReader();

	virtual bool readEntry();
	virtual int getNumFields() const { return _numFields; }
	virtual void getField(int numField, string &val) const;
	virtual void getField(int numField, CHRPOS &val) const; //this signaiture isn't const because it operates on an internal string for speed.
	virtual void getField(int fieldNum, char &val) const;
	virtual void appendField(int fieldNum, string &str) const;
	virtual const string &getHeader() const { return _header; }
	virtual bool hasHeader() const { return _fullHeaderFound; }
	virtual void setInHeader(bool val) { _inheader = val; }

protected:
	int _numFields;
	char _delimChar;
	string _header;
	bool _fullHeaderFound;
	int _currDataPos;
	string _sLine;
	int *_delimPositions;
	string _currChromStr;
	string _tempChrPosStr;
	int _lineNum;
	bool _inheader;
	bool detectAndHandleHeader();
	bool findDelimiters();

	// Resolve a field's [start, end) bounds from the compact position
	// table. Field offsets are kept in the reader's line-local coordinate
	// space, so all accessors normalize through these helpers rather
	// than indexing _delimPositions directly.
	int startOffset(int fieldNum) const;
	int endOffset(int fieldNum) const;
	int fieldLength(int fieldNum) const;

	//This is actually a very specialized function strictly for VCF
	//records to read and process extra information about Structural Variants.
	static const int VCF_TAG_FIELD = 7;
	CHRPOS getVcfSVlen();


};


#endif /* SINGLELINETEXTFILEREADER_H_ */
