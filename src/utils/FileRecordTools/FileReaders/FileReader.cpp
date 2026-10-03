#include <iostream>
#include <cstring>
#include "FileReader.h"
#include "BufferedStreamMgr.h"

FileReader::FileReader()
: _fileIdx(-1),
  _bufStreamMgr(NULL),
  _isFileOpen(false),
  _currChromId(-1)
{
}

FileReader::~FileReader() {
}

bool FileReader::open() {

	// The input stream is wired up by setInputStream() before open() is
	// called; if that hasn't happened we have no stream to read from.
	if (_isFileOpen) {
		return true;
	} else {
		fprintf(stderr, "Error: bad input stream.\n");
		exit(1);
	}
	return false; //can't reach this point, but it's here to satisfy the compiler.
}

void FileReader::close() {
}

bool FileReader::eof() const {
	return _bufStreamMgr->eof();
}
