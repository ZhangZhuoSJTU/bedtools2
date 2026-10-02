/*
 * ContextGroupBy.cpp
 *
 *  Created on: Mar 26, 2014
 *      Author: nek3d
 */


#include "ContextGroupBy.h"

ContextGroupBy::ContextGroupBy()
: _printFullCols(false),
  _ignoreCase(false)
{
	setSortedInput(true);
	_noEnforceCoordSort = true;
	setColumnOpsMethods(true);
	setExplicitBedOutput(true);

	//For columnOps, groupBy has default operation sum but no default column,
	///so we have to clear that.
	//
	_keyListOps->setColumns("");

}

ContextGroupBy::~ContextGroupBy()
{

}


bool ContextGroupBy::parseCmdArgs(int argc, char **argv, int skipFirstArgs)
{
	struct ArgHandler {
		const char *name;
		const char *alias;
		bool (ContextGroupBy::*handler)();
	};
	static const ArgHandler handlers[] = {
		{"-g", "-grp", &ContextGroupBy::handle_g},
		{"-inheader", NULL, &ContextGroupBy::handle_inheader},
		{"-outheader", NULL, &ContextGroupBy::handle_outheader},
		{"-header", NULL, &ContextGroupBy::handle_header},
		{"-full", NULL, &ContextGroupBy::handle_full},
		{"-ignorecase", NULL, &ContextGroupBy::handle_ignorecase},
	};
	const size_t numHandlers = sizeof(handlers) / sizeof(handlers[0]);

	for (_i=_skipFirstArgs; _i < argc; _i++) {
		if (isUsed(_i - _skipFirstArgs)) {
			continue;
		}

		for (size_t h = 0; h < numHandlers; h++) {
			if (strcmp(_argv[_i], handlers[h].name) == 0 ||
				(handlers[h].alias && strcmp(_argv[_i], handlers[h].alias) == 0)) {
				if (!(this->*handlers[h].handler)()) return false;
				break;
			}
		}
	}
	return ContextBase::parseCmdArgs(argc, argv, _skipFirstArgs);
}

bool ContextGroupBy::isValidState()
{
	// The user was required to have entered one or more columns
	if (_keyListOps->getColumns().empty()) {
		_errorMsg = "***** ERROR: -opCols parameter requires a value.";
		return false;
	}

	//default to stdin
	if (getNumInputFiles() == 0) {
		addInputFile("-");
	}
	//default grouping is cols 1,2,3
	if (_groupStr.empty()) _groupStr = "1,2,3";

	return ContextBase::isValidState();
}


bool ContextGroupBy::handle_g()
{
	if (_argc <= _i+1) {
		_errorMsg = "\n***** ERROR: -g option given, but columns to group not specified. *****";
		return false;
	}
	_groupStr = _argv[_i+1];
	markUsed(_i - _skipFirstArgs);
	_i++;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextGroupBy::handle_inheader()
{
	_inheader = true;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextGroupBy::handle_outheader() {
	setPrintHeader(true);
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextGroupBy::handle_header() {
	_inheader = true;
	setPrintHeader(true);
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextGroupBy::handle_full() {
	_printFullCols = true;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextGroupBy::handle_ignorecase() {
	_ignoreCase = true;
	markUsed(_i - _skipFirstArgs);
	return true;
}

const string &ContextGroupBy::getDefaultHeader() {
	//groupBy does not support multiple databases.
	FileRecordMgr *frm = _files[0];
	int numFields = frm->getNumFields();
	_defaultHeader.clear();
	ostringstream s;
	for (int i=1; i <= numFields; i++) {
		s << "col_";
		s << i;
		s << "\t";
	}
	_defaultHeader.append(s.str());
	//change last tab into newline
	_defaultHeader[_defaultHeader.size()-1] = '\n';
	return _defaultHeader;
}
