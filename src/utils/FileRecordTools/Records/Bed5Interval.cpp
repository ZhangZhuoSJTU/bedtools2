#include "Bed5Interval.h"
#include "SingleLineDelimTextFileReader.h"
#include <cstring>

Bed5Interval::Bed5Interval()
{

}

Bed5Interval::~Bed5Interval()
{

}

bool Bed5Interval::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	// A BED5 record extends the BED3 core with the name and score columns.
	bool parsedBase = Bed3Interval::initFromFile(fileReader);

	// Columns 3 and 4 (0-based) hold the name and score respectively.
	fileReader->getField(3, _name);
	fileReader->getField(4, _score);
	return parsedBase;
}

inline void print_record(const string& chrom_name, CHRPOS start, CHRPOS end, const string& name, const string& score, string& buf) {
	const char* buffer = buffer_printf("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS"\t%s\t%s", chrom_name.c_str(), start, end, name.c_str(), score.c_str());
	buf.append(buffer);
}

void Bed5Interval::print(string &outBuf) const
{
	print_record(_chrName, _startPos, _endPos, _name, _score, outBuf);
}

void Bed5Interval::print(string &outBuf, CHRPOS start, CHRPOS end) const
{
	print_record(_chrName, start, end, _name, _score, outBuf);
}

void Bed5Interval::print(string &outBuf, const string & start, const string & end) const
{
	Bed3Interval::print(outBuf, start, end);
	outBuf.append("\t");
	outBuf.append(_name);
	outBuf.append("\t");
	outBuf.append(_score);
}


void Bed5Interval::printNull(string &outBuf) const
{
	Bed3Interval::printNull(outBuf);
	outBuf.append("\t.\t-1", 5);
}

const string &Bed5Interval::getField(int fieldNum) const
{
	switch (fieldNum) {
	case NAME_COLUMN:
		return _name;
		break;
	case SCORE_COLUMN:
		return _score;
		break;

	default:
		return Bed3Interval::getField(fieldNum);
		break;
	}
}

bool Bed5Interval::isNumericField(int fieldNum) {
	switch (fieldNum) {
	case NAME_COLUMN:
		return false;
		break;
	case SCORE_COLUMN:
		return true;
		break;
	default:
		return Bed3Interval::isNumericField(fieldNum);
	}
}
