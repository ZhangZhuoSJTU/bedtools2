/*
 * RegressTestMemory.cpp
 *
 * Memory profiling, split out of RegressTest.cpp
 */

#include "RegressTest.h"
#include <cstdlib>
#include "SingleLineDelimTextFileReader.h"
#include "BufferedStreamMgr.h"
#include "lineFileUtilities.h"
#include "ParseTools.h"
#include <sys/stat.h>
#include <ctime>

bool RegressTest::startMemoryProfile(bool isBedops)
{
	//kick off ps as a background process to monitor memory usage
	string pidFilename = _tmpDirname + "pidFile.txt";
	string psCmd;
	if (!isBedops) {
		psCmd = "while [ 1 ] ; do ps au | grep bedtools | grep ";
	} else {
		psCmd = "while [ 1 ] ; do ps au | grep bedops | grep ";
	}
	psCmd += _userName;
	psCmd += " | grep -v grep | grep -v regresstest | grep -v \"sh -c\" >> ";
	psCmd += _memoryLogfilename;
	psCmd += " ; sleep .5 ; done & echo $! > ";
	psCmd += _tmpDirname + "pidFile.txt";

	if (system(psCmd.c_str()) != 0) {
		fprintf(stderr, "Error: unable to launch top as background process for memory profiling.\n");
		return false;
	}

	// open and parse the pidFile to get the process id of the memory profiling background process.
	FILE *fp = fopen(pidFilename.c_str(), "r");
	char sLine[4192];
	memset(sLine, 0, 4192);
	fgets(sLine, 4192, fp);
	fclose(fp);

	_currPidOfMemoryLogging = sLine;
	return true;
}

bool RegressTest::endMemoryProfile()
{

	string killCmd = "kill -9 ";
	killCmd += _currPidOfMemoryLogging;
	if (system(killCmd.c_str()) != 0) {
		fprintf(stderr, "Error: failed to kill process id %s\n", _currPidOfMemoryLogging.c_str());
		return false;
	}

	return true;
}

bool RegressTest::calcMemoryStats()
{
	//read and parse the memory log file, calc basic stats: max, mean, median.
	FILE *fp = fopen(_memoryLogfilename.c_str(), "r");
	if (fp == NULL) {
		fprintf(stderr, "Error: Unable to read memory profile log %s.\n", _memoryLogfilename.c_str());
		return false;
	}

	char sLine[4192];
	int numValidLines=0;
	char field1[2048];
	char field2[2048];
	char field3[2048];
	char field4[2048];
	char field5[2048];

	int totalMem= 0;
	int currMem = 0;
	int maxMem = 0;
	while(!feof(fp)) {
		memset(sLine, 0, 4192);
		memset(field1, 0, 2048);
		memset(field2, 0, 2048);
		memset(field3, 0, 2048);
		memset(field4, 0, 2048);
		memset(field5, 0, 2048);

		currMem = 0;
		fgets(sLine, 4192, fp);
		int len = strlen(sLine);
		bool isValidLine = false;
		for (int i=0; i < len; i++) {
			if (!isspace(sLine[i])) {
				isValidLine = true;
				break;
			}
		}
		if (isValidLine) {
			numValidLines++;
			sscanf(sLine, "%s %s %s %s %s", field1, field2, field3, field4, field5);

			//now field 5 has the number of kilobytes the process used at that moment.
			currMem = atoi(field5);
			if (currMem == 0) {
				//atoi failed. something wrong with input
				fprintf(stderr, "Error: bad field where memory usage expected: %s.\n", field5);
				fclose(fp);
				return false;
			}
			maxMem = max(currMem, maxMem);
			totalMem += currMem;
		}
	}

	int avgMem = numValidLines > 0 ? totalMem / numValidLines : -1;

	fclose(fp);
	if (avgMem > 0) {
		fprintf(_fpReportFile, "\tMemory uage: max = %dkb\tavg = %dkb\n", maxMem, avgMem);
	} else {
		fprintf(_fpReportFile, "\tMemory usage: TOO FAST TO PROFILE!\n");
	}

	//cleanup: erase memoryLogfile, as we wish to append to a blank file when next used.
	if (remove(_memoryLogfilename.c_str()) != 0) {
		fprintf(stderr, "Error: couldn't delete old memory logfile %s.\n", _memoryLogfilename.c_str());
		return false;
	}

	return true;
}
