/*
 * ContextBaseHandlers.cpp
 *
 *  Command line option handlers split out of ContextBase.cpp
 */

#include "ContextBase.h"
#include "Random.h"
#include <unistd.h>
#include <sys/types.h>
#include <cctype>

bool ContextBase::handle_bed()
{
	setExplicitBedOutput(true);
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_fbam()
{
	setUseFullBamTags(true);
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_g()
{
	if (_argc <= _i+1) {
		_errorMsg = "\n***** ERROR: -g option given, but no genome file specified. *****";
		return false;
	}
	openGenomeFile(_argv[_i+1]);
	markUsed(_i - _skipFirstArgs);
	_i++;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_h()
{
	setShowHelp(true);
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_header()
{
	setPrintHeader(true);
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_i()
{
	if (_argc <= _i+1) {
		_errorMsg = "\n***** ERROR: -i option given, but no input file specified. *****";
		return false;
	}
	addInputFile(_argv[_i+1]);
	markUsed(_i - _skipFirstArgs);
	_i++;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_n()
{
	if (_argc <= _i+1) {
		_errorMsg = "\n***** ERROR: -n option given, but no number of output records specified. *****";
		return false;
	}
	setNumOutputRecords(atoi(_argv[_i + 1]));
	markUsed(_i - _skipFirstArgs);
	_i++;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_nobuf()
{
	setUseBufferedOutput(false);
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_iobuf()
{
	if (_argc <= _i+1) {
		_errorMsg = "\n***** ERROR: -iobuf option given, but size of input buffer not specified. *****";
		return false;
	}
	if (!parseIoBufSize(_argv[_i + 1])) return false;
	markUsed(_i - _skipFirstArgs);
	_i++;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_seed()
{
	if (_argc <= _i+1) {
		_errorMsg = "\n***** ERROR: -seed option given, but no seed specified. *****";
		return false;
	}
	_hasConstantSeed = true;
	_seed  = atoi(_argv[_i+1]);
	rand_set_seed(_seed);
	markUsed(_i - _skipFirstArgs);
	_i++;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_split()
{
    setObeySplits(true);
     markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_sorted()
{
	setSortedInput(true);
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_ubam()
{
    setUncompressedBam(true);
    markUsed(_i - _skipFirstArgs);
	return true;
}


// Methods specific to column operations.
// for col ops, -c is the string of columns upon which to operate
bool ContextBase::handle_c()
{
	if (!hasColumnOpsMethods()) {
		return false;
	}
    if ((_i+1) < _argc) {
        _keyListOps->setColumns(_argv[_i + 1]);
        markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
        return true;
    }
    return false;
}


// for col ops, -o is the string of operations to apply to the columns (-c)
bool ContextBase::handle_o()
{
	if (!hasColumnOpsMethods()) {
		return false;
	}
    if ((_i+1) < _argc) {
    	 _keyListOps->setOperations(_argv[_i + 1]);
        markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
        return true;
    }
    return false;
}

bool ContextBase::handle_prec()
{
	if (!hasColumnOpsMethods()) {
		return false;
	}
    if ((_i+1) < _argc) {
    	int prec = atoi(_argv[_i + 1]);
    	if (prec < 1) {
    		_errorMsg += "\n***** ERROR: -prec must be followed by a positive integer. Exiting. *****";
    		return false;
    	}
    	 _reportPrecision = prec;
        markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
        return true;
    }
	_errorMsg += "\n***** ERROR: -prec must be followed by a positive integer. Exiting. *****";
    return false;
}



// for col ops, -null is a NULL value assigned
// when no overlaps are detected.
bool ContextBase::handle_null()
{
	if (!hasColumnOpsMethods()) {
		return false;
	}
    if ((_i+1) < _argc) {
    	 _keyListOps->setNullValue(_argv[_i + 1]);
        markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
        return true;
    }
    return false;
}

//for col ops, delimStr will appear between each item in
//a collapsed but delimited list.
bool ContextBase::handle_delim()
{
	if (!hasColumnOpsMethods()) {
		_errorMsg = "\n***** ERROR: Can't set delimiter for tools without column operations. Exiting. *****";
		return false;
	}
    if ((_i+1) < _argc) {
    	 _keyListOps->setDelimStr(_argv[_i + 1]);
        markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
    }
    return true;
}

bool ContextBase::handle_sortout()
{
	setSortOutput(true);
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_nonamecheck()
{
	setNameCheckDisabled(true);
	markUsed(_i - _skipFirstArgs);
	return true;
}
