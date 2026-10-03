#include "PlusFields.h"
#include "SingleLineDelimTextFileReader.h"

PlusFields::PlusFields()
{
}

bool PlusFields::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	// The "plus" fields are every column beyond the fixed record prefix
	// (_numOffsetFields). Size our storage to hold exactly those trailing
	// columns, growing/shrinking only when the count changes.
	size_t numPlusFields = fileReader->getNumFields() - _numOffsetFields;

	if (size() != numPlusFields) {
		_fields.resize(numPlusFields);
	}

	// Copy each trailing column, translating our local index into the
	// reader's absolute 0-based column index.
	for (size_t i = 0; i < numPlusFields; i++) {
		fileReader->getField((int)(i + _numOffsetFields), _fields[i]);
	}
	return true;
}

void PlusFields::clear() {
	//don't destroy the strings if we don't have to. Just clear their memory.
	for (int i=0; i < (int)_fields.size(); i++) {
		_fields[i].clear();
	}
}

const string &PlusFields::getField(int fieldNum) const
{
	return _fields[fieldNum - _numOffsetFields - 1];
}

void PlusFields::printFields(string &outBuf) const {
	for (size_t i=0; i < size(); i++) {
		outBuf.append(_fields[i]);
		if (i < size() - 1) outBuf.append("\t");
	}
}
