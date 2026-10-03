/*
 * RecordOutputMgr.h
 *
 *  Created on: May 28, 2013
 *      Author: nek3d
 */

#ifndef RECORDOUTPUTMGR_H_
#define RECORDOUTPUTMGR_H_

#include "ContextBase.h"
#include "RecordKeyVector.h"
#include "api/BamWriter.h"

using namespace std;

class BlockMgr;

// Central sink for all tool output. Buffers textual (BED/VCF/GFF) output and
// flushes it to stdout in blocks, or writes native BAM through a BamWriter.
// Knows how to render the many per-tool/per-option output layouts (intersect,
// closest, map, merge, groupBy, sample, ...).
class RecordOutputMgr {
public:
	RecordOutputMgr();
	~RecordOutputMgr();

	//The init method must be called after all the input files are open.
	void init(ContextBase *context);

	// Emit a single record (optionally paired with a trailing value/column).
	void printRecord(Record *record);
	void printRecord(RecordKeyVector &keyList);
	void printRecord(Record *record, const string & value);

	// Emit any pending file header exactly once, before the first record.
	void checkForHeader();

	// Emit a query record together with its closest hit(s) and, optionally,
	// the signed/absolute distances to each hit.
	void printClosest(RecordKeyVector &keyList, const vector<CHRPOS> *dists = NULL);

	// Append the field (tab) and line (newline) separators to the output buffer.
	void tab() { _outBuf.append(FIELD_SEPARATOR); }
	void newline() { _outBuf.append(RECORD_SEPARATOR); }

private:
	typedef enum { NOT_BAM, BAM_AS_BAM, BAM_AS_BED} printBamType;

	ContextBase *_context;
	bool _printable;
	BamTools::BamWriter *_bamWriter;
	RecordKeyVector *_currBamBlockList;

	string _outBuf;

	// Delimiters used when assembling a line of textual (non-BAM) output.
	static constexpr const char *FIELD_SEPARATOR = "\t";
	static constexpr const char *RECORD_SEPARATOR = "\n";

	BlockMgr *_bamBlockMgr;
	string _afterVal; //to store values to be printed after record, such as column operations.
	//some helper functions to neaten the code.

	// Emit the null/blank form of the query or database record type.
	void null(bool queryType, bool dbType);

	// Core record emitter shared by the public printRecord overloads.
	void printRecord(RecordKeyVector &keyList, RecordKeyVector *blockList);

	// Print just a key record, with or without overriding start/end coordinates.
	void printKey(const Record *key);
	void printKey(const Record *key, const string & start, const string & end);
	void printKey(const Record *key, CHRPOS start, CHRPOS end);

	// Emit the leading database-file-id column when multiple databases are used.
	void addDbFileId(int fileId);

	// Print the key record and report whether it was emitted as native BAM.
	bool printKeyAndTerminate(RecordKeyVector &keyList);

	// Route a possibly-BAM record to the BAM writer or the text buffer.
	printBamType printBamRecord(RecordKeyVector &keyList, bool bamOutputOnly = false);

	// Emit one query/hit pair (per -wa/-wb/-wo/-loj) or a per-query summary.
	void reportOverlapDetail(const Record *keyRecord, const Record *hitRecord, int hitIdx = 0);
	void reportOverlapSummary(RecordKeyVector &keyList);

	static const unsigned int MAX_OUTBUF_SIZE = 16384; //16 K

	// When buffered output is in use, the buffer is flushed once it reaches this
	// fraction of MAX_OUTBUF_SIZE (i.e. 90% full).
	static constexpr double FLUSH_FILL_FRACTION = 0.9;

	// If we are using buffered output, only flush the output buffer if it's least
	// 90% full. If we're not using buffered output, flush if it's not empty
	bool needsFlush() const {
		return ((_context->getUseBufferedOutput() &&_outBuf.size() >= MAX_OUTBUF_SIZE * FLUSH_FILL_FRACTION) ||
				(!_context->getUseBufferedOutput() && !_outBuf.empty()));
	}
	void flush();
};

#endif /* RECORDOUTPUTMGR_H_ */
