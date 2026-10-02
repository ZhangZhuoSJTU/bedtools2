/*
 * ContextMerge.cpp
 *
 *  Created on: Mar 26, 2014
 *      Author: nek3d
 */


#include "ContextMerge.h"

ContextMerge::ContextMerge()
{
	setSortedInput(true);
	setUseMergedIntervals(true);
	setColumnOpsMethods(true);
	setExplicitBedOutput(true);

	//merge has no default columnOps the way map does, so we'll need to clear those.
	_keyListOps->setColumns("");
	_keyListOps->setOperations("");

}

ContextMerge::~ContextMerge()
{

}


bool ContextMerge::parseCmdArgs(int argc, char **argv, int skipFirstArgs)
{
	struct ArgHandler {
		const char *name;
		const char *alias;
		bool (ContextMerge::*handler)();
	};
	static const ArgHandler handlers[] = {
		{"-n", NULL, &ContextMerge::handle_n},
		{"-nms", NULL, &ContextMerge::handle_nms},
		{"-scores", NULL, &ContextMerge::handle_scores},
		{"-delim", NULL, &ContextMerge::handle_delim},
		{"-d", NULL, &ContextMerge::handle_d},
		{"-s", NULL, &ContextMerge::handle_s},
		{"-S", NULL, &ContextMerge::handle_S},
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

bool ContextMerge::isValidState()
{
	// Special: The merge program does not have default
	//column operations, so if none were entered, disable column ops.
	if (_keyListOps->getColumns().empty() && _keyListOps->getOperations().empty()) {
		setColumnOpsMethods(false);
		delete _keyListOps;
		_keyListOps = NULL;
	}

	//default to stdin
	if (getNumInputFiles() == 0) 
	{
		if (!isatty(STDIN_FILENO))
		{
			addInputFile("-");
		}
	}
	if (!ContextBase::isValidState()) {
		return false;
	}

	//
	// Tests for stranded merge
	//
	if (_desiredStrand != FileRecordMergeMgr::ANY_STRAND) { // requested stranded merge
		return strandedToolSupported();
	}

	return true;
}


bool ContextMerge::handle_d() {
    if ((_i+1) < _argc) {
    	if (isNumeric(_argv[_i+1])) {
			int dist = (int)str2chrPos(_argv[_i+1]);
			
			_maxDistance = dist;
	    	markUsed(_i - _skipFirstArgs);
	        _i++;
	        markUsed(_i - _skipFirstArgs);
			return true;
    	}
    }
	_errorMsg = "\n***** ERROR: -d option must be followed by an integer value *****";
	return false;
}

bool ContextMerge::handle_n()
{
	markUsed(_i - _skipFirstArgs);
	_errorMsg = "\n***** ERROR: -n option is deprecated. Please see the documentation for the -c and -o column operation options. *****";
    return false;
}

bool ContextMerge::handle_nms()
{
	markUsed(_i - _skipFirstArgs);
	_errorMsg = "\n***** ERROR: -nms option is deprecated. Please see the documentation for the -c and -o column operation options. *****";
    return false;

}


bool ContextMerge::handle_scores()
{
	// No longer supporting this deprecated option.
	markUsed(_i - _skipFirstArgs);
	_errorMsg = "\n***** ERROR: -scores option is deprecated. Please see the documentation for the -c and -o column operation options. *****";
    return false;
}

bool ContextMerge::handle_s() {
	_desiredStrand = FileRecordMergeMgr::SAME_STRAND_EITHER;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextMerge::handle_S() {
    if ((_i+1) < _argc) {
    	bool validChar = false;
    	if (_argv[_i+1][0] == '+') {
			_desiredStrand = FileRecordMergeMgr::SAME_STRAND_FORWARD;
			validChar = true;
    	} else if (_argv[_i+1][0] == '-') {
    		validChar = true;
			_desiredStrand = FileRecordMergeMgr::SAME_STRAND_REVERSE;
    	}
    	if (validChar) {
			markUsed(_i - _skipFirstArgs);
			_i++;
			markUsed(_i - _skipFirstArgs);
			return true;
    	}
    }
	_errorMsg = "\n***** ERROR: -S option must be followed by + or -. *****";
	return false;
}
