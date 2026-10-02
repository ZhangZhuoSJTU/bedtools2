/*
 * RecDistList.cpp
 *
 * RecDistList methods, split out of CloseSweep.cpp
 */

#include "CloseSweep.h"

RecDistList::RecDistList(int maxSize)
:  _kVal(maxSize),
   _empty(true),
   _currNumIdxs(0),
   _totalRecs(0)

    {

	_allRecs.resize(_kVal);
	for (int i=0; i < _kVal; i++) {
		_allRecs[i] = new elemsType();
	}
	_distIndex = new indexType[_kVal];

	clear();
}

RecDistList::~RecDistList() {
	for (int i=0; i < _kVal; i++) {
		delete _allRecs[i];
	}
	delete _distIndex;
	_distIndex = NULL;

}

void RecDistList::clear() {
	for (int i=0; i < _kVal; i++) {
		_allRecs[i]->clear();
		_allRecs[i]->reserve(16);
	}
	for (int i=0; i < _kVal; i++) {
		_distIndex[i].first = -1;
		_distIndex[i].second = -1;
	}
	_currNumIdxs = 0;
	_empty = true;
	_totalRecs = 0;
}

bool RecDistList::addRec(CHRPOS dist, Record *record, chromDirType chromDir) {
	CHRPOS newPos = 0;
	bool mustAppend = false;
	int useElemIdx = 0;

	if (dist > getMaxDist()) {
		//dist is bigger than any currently contained.
		if (_currNumIdxs == _kVal) {
			return false;
		}
		mustAppend = true;
		useElemIdx = _currNumIdxs;
		newPos = _currNumIdxs;
	}
	if (find(dist, newPos)) {
		_allRecs[_distIndex[newPos].second]->push_back(elemPairType(chromDir, record));
		_totalRecs++;
		return true;
	}
	if (!mustAppend) {
		//smaller than maxDist, and doesn't currently exist. must insert
		//newPos now is the insertion point.
		CHRPOS startShiftPos = 0;
		if (_currNumIdxs == _kVal) {
			//already full. must remove oldest max.
			//determine which vector it was using
			//so we can re-use it.
			startShiftPos = _kVal-1;
			useElemIdx = _distIndex[startShiftPos].second;
		} else {
			//can add a new element
			startShiftPos = _currNumIdxs;
			useElemIdx = _currNumIdxs;
			_currNumIdxs++;
		}
		for (CHRPOS i=startShiftPos; i > newPos; i--) {
			_distIndex[i].first = _distIndex[i-1].first;
			_distIndex[i].second = _distIndex[i-1].second;
		}
	} else {
		_currNumIdxs++;
	}
	_allRecs[useElemIdx]->clear();
	_allRecs[useElemIdx]->reserve(16);
	_allRecs[useElemIdx]->push_back(elemPairType(chromDir, record));

	_distIndex[newPos].first = (int)dist;
	_distIndex[newPos].second = useElemIdx;
	_empty = false;
	_totalRecs++;
	return true;
}

//if true, pos will be the idx the distance is at.
//if false, pos will be the idx to insert at.
bool RecDistList::find(CHRPOS dist, CHRPOS &pos) const {
	CHRPOS lbound=0, ubound=_currNumIdxs-1, currVal =0;
	pos = 0;
	while(lbound <= ubound)
	{
		pos = (lbound + ubound) / 2;
		currVal = _distIndex[pos].first;
		if (currVal == dist) {
			return true;
		}
		if (dist > currVal) {
			lbound = pos + 1;
		} else {
			ubound = pos -1;
		}
	}
	pos = (ubound == -1 ) ? 0 : (lbound == _currNumIdxs ? _currNumIdxs : lbound);
	return false;
}

CHRPOS RecDistList::getMaxLeftEndPos() const {

	if (_empty) return -1;

	CHRPOS maxDist =_distIndex[_currNumIdxs-1].first;

  const elemsType *elems = _allRecs[_distIndex[_currNumIdxs-1].second];
	for (int i=0; i < (int)elems->size(); i++) {
    const elemPairType & elem = (*elems)[i];
		if (elem.first == LEFT) {
			return maxDist;
		}
	}
	return -1;

}
