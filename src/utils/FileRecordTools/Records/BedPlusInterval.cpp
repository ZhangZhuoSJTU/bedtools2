#include "BedPlusInterval.h"
#include "SingleLineDelimTextFileReader.h"

BedPlusInterval::BedPlusInterval()
:  _numFixedFields(defaultNumFixedFields),
  _numPrintFields(0)
{
	_plusFields.setNumOffsetFields(defaultNumFixedFields);
}

void BedPlusInterval::setNumFixedFields(int numFields) {
	_numFixedFields = numFields;
	_plusFields.setNumOffsetFields(numFields);
}


bool BedPlusInterval::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	// Always parse the BED3 core first; the "plus" columns are layered on
	// top of whatever fixed-field width this record was configured with.
	bool parsedBase = Bed3Interval::initFromFile(fileReader);

	// When the fixed portion is wider than BED3, pull the BED6 columns
	// (name/score/strand at 0-based indices 3-5) before the trailing
	// variable-width plus fields.
	if (_numFixedFields != defaultNumFixedFields) {
		fileReader->getField(3, _name);
		fileReader->getField(4, _score);
		fileReader->getField(5, _strand);
		adjustStrandVal();
	}
	_plusFields.initFromFile(fileReader);
	return parsedBase;

}


void BedPlusInterval::clear() {
	Bed3Interval::clear();
	_plusFields.clear();
}

void BedPlusInterval::print(string &outBuf) const
{
	Bed3Interval::print(outBuf);
	outBuf.append("\t");
	printBed6PlusFields(outBuf);
	_plusFields.printFields(outBuf);
}

void BedPlusInterval::print(string &outBuf, CHRPOS start, CHRPOS end) const
{
	Bed3Interval::print(outBuf, start, end);
	outBuf.append("\t");
	printBed6PlusFields(outBuf);
	_plusFields.printFields(outBuf);
}

void BedPlusInterval::print(string &outBuf, const string & start, const string & end) const
{
	Bed3Interval::print(outBuf, start, end);
	outBuf.append("\t");
	printBed6PlusFields(outBuf);
	_plusFields.printFields(outBuf);
}


void BedPlusInterval::printNull(string &outBuf) const
{
	Bed3Interval::printNull(outBuf);
	printBed6PlusNullFields(outBuf);
	for (int i=_numFixedFields; i < _numPrintFields; i++) {
		outBuf.append("\t.");
	}
}

const string &BedPlusInterval::getField(int fieldNum) const
{
	if (fieldNum > _numFixedFields) {
		return _plusFields.getField(fieldNum);
	} else if (fieldNum == 4 && _numFixedFields >=4) {
		return _name;
	} else if (fieldNum == 5 && _numFixedFields >=5) {
		return _score;
	}
	else if (fieldNum == 6 && _numFixedFields >=6) {
		return _strand;
	}
	return Bed3Interval::getField(fieldNum);
}

bool BedPlusInterval::isNumericField(int fieldNum) {

	//
	// TBD: There is no currently no good way to guarantee / enforce whether
	// fields after the 3rd are numeric, so for now we'll give the user the
	// benefit of the doubt on those.
	//
	if (fieldNum > defaultNumFixedFields) {
		return true;
	}
	return Bed3Interval::isNumericField(fieldNum);
}

void BedPlusInterval::printBed6PlusFields(string &outBuf) const {
	if (_numFixedFields != defaultNumFixedFields) {
		outBuf.append(_name);
		outBuf.append("\t");
		outBuf.append(_score);
		outBuf.append("\t");
		outBuf.append(_strand);
		outBuf.append("\t");
	}
}

void BedPlusInterval::printBed6PlusNullFields(string &outBuf) const {
	if (_numFixedFields != defaultNumFixedFields) {
		outBuf.append("\t.\t.\t.");
	}

}
