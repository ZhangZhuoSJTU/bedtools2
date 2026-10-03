/*
 * Bed6Interval.h
 *
 *  Created on: Nov 13, 2012
 *      Author: nek3d
 */

#ifndef BED6INTERVAL_H_
#define BED6INTERVAL_H_


#include "Bed3Interval.h"

class SingleLineDelimTextFileReader;

class Bed6Interval : public Bed3Interval {
public:
	friend class FreeList<Bed6Interval>;

	Bed6Interval();
	virtual bool initFromFile(SingleLineDelimTextFileReader *);
	virtual void print(string &outBuf) const;
	virtual void print(string &outBuf, CHRPOS start, CHRPOS end) const;
	virtual void print(string &outBuf, const string & start, const string & end) const;
	virtual void printNull(string &outBuf) const;
	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::BED6_RECORD_TYPE; }

	virtual const string &getField(int fieldNum) const;
	virtual int getNumFields() const  { return 6; }
	static bool isNumericField(int fieldNum);

	// 1-based column positions addressed by the generic field accessors
	// (getField / isNumericField), distinct from the 0-based indices used
	// when pulling raw fields out of the reader in initFromFile.
	enum {
		NAME_COLUMN   = 4,
		SCORE_COLUMN  = 5,
		STRAND_COLUMN = 6
	};


protected:
	virtual ~Bed6Interval();
};

#endif /* BED6INTERVAL_H_ */
