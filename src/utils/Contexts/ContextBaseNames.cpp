/*
 * ContextBaseNames.cpp
 *
 * Chromosome naming convention checks and warnings, split out of ContextBase.cpp
 */

#include "ContextBase.h"
#include <cctype>

void ContextBase::testNameConventions(const Record *record) {
	//Do nothing if using the -nonamecheck option,
	//warning already given, or record is unmapped BAM record
	if (getNameCheckDisabled() || _nameConventionWarningTripped || record->isUnmapped()) return;

	int fileIdx = record->getFileIdx();

	//
	// First test whether chr in chrom names match
	//

	bool hasChr = record->hasChrInChromName();
	testType testChrVal = fileHasChrInChromNames(fileIdx);

	if (testChrVal == UNTESTED) {
		_fileHasChrInChromNames[fileIdx] = hasChr ? YES : NO;
	}
	if ((_allFilesHaveChrInChromNames == YES && !hasChr) || (_allFilesHaveChrInChromNames == NO && hasChr)) {
		nameConventionWarning(record, _fileNames[fileIdx], " has inconsistent naming convention for record:\n");
	}

	if (_allFilesHaveChrInChromNames == UNTESTED) {
		_allFilesHaveChrInChromNames = hasChr ? YES : NO;
	}


	//
	// Now test whether leading zero in chrom names match
	//


	bool zeroVal = record->hasLeadingZeroInChromName(hasChr);
	testChrVal = fileHasLeadingZeroInChromNames(fileIdx);
	if (testChrVal == UNTESTED) {
		_fileHasLeadingZeroInChromNames[fileIdx] = zeroVal ? YES : NO;
	}
	if ((_allFileHaveLeadingZeroInChromNames == YES && !zeroVal) || (_allFileHaveLeadingZeroInChromNames == NO && zeroVal)) {
		nameConventionWarning(record, _fileNames[fileIdx], " has a record where naming convention (leading zero) is inconsistent with other files:\n");
	}

	if (_allFileHaveLeadingZeroInChromNames == UNTESTED) {
		_allFileHaveLeadingZeroInChromNames = zeroVal ? YES : NO;
	}
}

ContextBase::testType ContextBase::fileHasChrInChromNames(int fileIdx) {
	conventionType::iterator iter = _fileHasChrInChromNames.find(fileIdx);
	if (iter == _fileHasChrInChromNames.end()) {
		return UNTESTED;
	}
	return iter->second;
}

ContextBase::testType ContextBase::fileHasLeadingZeroInChromNames(int fileIdx) {
	conventionType::iterator iter = _fileHasLeadingZeroInChromNames.find(fileIdx);
	if (iter == _fileHasLeadingZeroInChromNames.end()) {
		return UNTESTED;
	}
	return iter->second;
}



void ContextBase::warn(const Record *record, const string str1, const string str2, const string str3) {
	string msg;
	setErrorMsg(msg, true, record, str1, str2, str3);
	cerr << msg << endl;
}

void ContextBase::die(const Record *record, const string str1, const string str2, const string str3) {
	string msg;
	setErrorMsg(msg, false, record, str1, str2, str3);
	cerr << msg << endl;
	exit(1);
}

void ContextBase::setErrorMsg(string &msg, bool onlyWarn, const Record * record, string str1, const string str2, const string str3) {
	if (onlyWarn) {
		msg = "\n***** WARNING: ";
	} else {
		msg = "\n***** ERROR: ";
	}
	msg.append(str1);
	msg.append(str2);
	msg.append(str3);
	msg.append(" Exiting...\n");
	if (record != NULL) {
		record->print(msg);
	}
}

void ContextBase::nameConventionWarning(const Record *record, const string &filename, const string &message)
{
	 _nameConventionWarningMsg = "***** WARNING: File ";
	 _nameConventionWarningMsg.append(filename);
	 _nameConventionWarningMsg.append(message);
	 record->print(_nameConventionWarningMsg);
	 _nameConventionWarningMsg.append("\n");
	 _nameConventionWarningTripped = true;

	 cerr << _nameConventionWarningMsg << endl;
}
