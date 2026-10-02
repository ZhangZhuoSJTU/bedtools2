#include "BedtoolsDriver.h"
#include <iostream>
#include "RecordOutputMgr.h"
//contexts
#include "ContextIntersect.h"
#include "ContextMerge.h"
#include "ContextMap.h"
#include "ContextSample.h"
#include "ContextJaccard.h"
#include "ContextClosest.h"
#include "ContextComplement.h"
#include "ContextSubtract.h"
#include "ContextSpacing.h"
#include "ContextCoverage.h"
#include "ContextSummary.h"

//tools
#include "intersectFile.h"
#include "mapFile.h"
#include "closestFile.h"
#include "mergeFile.h"
#include "jaccard.h"
#include "subtractFile.h"
#include "sampleFile.h"
#include "spacingFile.h"
#include "fisher.h"
#include "coverageFile.h"
#include "complementFile.h"
#include "groupBy.h"
#include "summaryFile.h"

template <class C>
static ContextBase *newContext() { return new C(); }

template <class C, class T>
static ToolBase *newTool(ContextBase *context) { return new T(static_cast<C *>(context)); }

struct ToolEntry {
	const char *name;
	ContextBase *(*makeContext)();
	ToolBase *(*makeTool)(ContextBase *);
};

static const ToolEntry toolTable[] = {
	{"intersect", newContext<ContextIntersect>, newTool<ContextIntersect, IntersectFile>},
	{"map", newContext<ContextMap>, newTool<ContextMap, MapFile>},
	{"closest", newContext<ContextClosest>, newTool<ContextClosest, ClosestFile>},
	{"merge", newContext<ContextMerge>, newTool<ContextMerge, MergeFile>},
	{"jaccard", newContext<ContextJaccard>, newTool<ContextJaccard, Jaccard>},
	{"subtract", newContext<ContextSubtract>, newTool<ContextSubtract, SubtractFile>},
	{"sample", newContext<ContextSample>, newTool<ContextSample, SampleFile>},
	{"spacing", newContext<ContextSpacing>, newTool<ContextSpacing, SpacingFile>},
	{"fisher", newContext<ContextFisher>, newTool<ContextFisher, Fisher>},
	{"coverage", newContext<ContextCoverage>, newTool<ContextCoverage, CoverageFile>},
	{"complement", newContext<ContextComplement>, newTool<ContextComplement, ComplementFile>},
	{"groupby", newContext<ContextGroupBy>, newTool<ContextGroupBy, GroupBy>},
	{"summary", newContext<ContextSummary>, newTool<ContextSummary, SummaryFile>},
};

static const ToolEntry *findTool(const string &name) {
	for (size_t i = 0; i < sizeof(toolTable) / sizeof(toolTable[0]); i++) {
		if (name == toolTable[i].name) return &toolTable[i];
	}
	cerr << "Error: Tool " << name << " is not supported. Exiting..." << endl;
	exit(1);
}

BedtoolsDriver::BedtoolsDriver()
: _hadError(false) {
	for (size_t i = 0; i < sizeof(toolTable) / sizeof(toolTable[0]); i++) {
		_supported.insert(toolTable[i].name);
	}
}


bool BedtoolsDriver::supports(const string &tool) {
	supportType::iterator iter = _supported.find(tool);
	return (iter != _supported.end());
}

bool BedtoolsDriver::subMain(int argc, char **argv) 
{
	_subCmd = argv[1];
	ContextBase *context = getContext();

	//process all command line arguments, check for valid usage,
	//show help and error messages if needed.
	if (!context->testCmdArgs(argc - 1, argv + 1)) {
		_hadError = context->errorEncountered();
		_errors = context->getErrorMessages();
		delete context;
		return false;
	}

	//establish which tool we're using (intersect, map, closest, etc).
	//initialize it.
	ToolBase *tool = getTool(context);
	if (!tool->init()) {
		delete context;
		return false;
	}

	//set-up output manager.
	RecordOutputMgr *outputMgr = new RecordOutputMgr();
	outputMgr->init(context);

	//process input
	RecordKeyVector hits;
	while (tool->findNext(hits)) {
		tool->processHits(outputMgr, hits);
		tool->cleanupHits(hits);
	}
	tool->finalizeCalculations();
	tool->giveFinalReport(outputMgr);
	delete outputMgr;
	delete tool;
	delete context;
	return true;
}


ContextBase *BedtoolsDriver::getContext()
{
	return findTool(_subCmd)->makeContext();
}

ToolBase *BedtoolsDriver::getTool(ContextBase *context)
{
	return findTool(_subCmd)->makeTool(context);
}
