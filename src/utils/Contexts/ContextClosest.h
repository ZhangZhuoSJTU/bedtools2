/*
 * ContextClosest.h
 *
 *  Created on: Sep 25, 2014
 *      Author: nek3d
 */


#ifndef CONTEXTCLOSEST_H_
#define CONTEXTCLOSEST_H_

#include "ContextIntersect.h"

// Context for the "closest" tool. It reuses all of the intersect options and
// adds options that control distance reporting (-d/-D), which neighbors to
// consider (-io/-iu/-id/-fu/-fd), tie handling (-t), multi-database behavior
// (-mdb), and how many closest hits to report (-k).
class ContextClosest : public ContextIntersect {
public:
	ContextClosest();
	virtual ~ContextClosest();
	virtual bool parseCmdArgs(int argc, char **argv, int skipFirstArgs);
    virtual bool hasIntersectMethods() const { return true; }
    virtual bool isValidState();

    // Default number of closest hits reported per query when -k is not given.
    static const int DEFAULT_NUM_CLOSEST_HITS = 1;

    bool hasTieMode() const { return _haveTieMode; }           // -t was given
    bool ignoreOverlaps() const { return _ignoreOverlaps; }    // -io
    bool ignoreUpstream() const { return _ignoreUpstream; }    // -iu
    bool ignoreDownstream() const { return _ignoreDownstream; }// -id
    bool forceUpstream() const { return _forceUpstream; }      // -fu
    bool forceDownstream() const { return _forceDownstream; }  // -fd
    bool reportDistance() const { return _reportDistance; }    // -d or -D
    bool signDistance() const { return _signDistance; }        // -D (signed)
    bool hasStrandedDistMode() const { return _haveStrandedDistMode; } // -D
    bool diffNames() const { return _diffNames; }              // -N
    int getNumClosestHitsWanted() const { return _numClosestHitsWanted; } // -k

    // Which end of a tie to keep when several hits are equidistant (-t).
    typedef enum { FIRST_TIE, LAST_TIE, ALL_TIES} tieModeType;
    tieModeType getTieMode() const { return _tieMode; }

    // Reference used to sign distances under -D (query "a", database "b", or
    // the reference-relative default).
    typedef enum { REF_DIST, A_DIST, B_DIST} strandedDistanceModeType;
    strandedDistanceModeType getStrandedDistMode() const { return _strandedDistMode; }

    // Whether multiple databases are searched independently or as one (-mdb).
    typedef enum { EACH_DB, ALL_DBS } multiDbModeType;
    multiDbModeType getMultiDbMode() const { return _multiDbMode; }

private:
    bool _haveTieMode;
    bool _ignoreOverlaps;
    bool _ignoreUpstream;
    bool _ignoreDownstream;
    bool _forceUpstream;
    bool _forceDownstream;
    bool _reportDistance;
    bool _signDistance;
    bool _haveStrandedDistMode;
    bool _diffNames;
    tieModeType _tieMode;
    strandedDistanceModeType _strandedDistMode;
    multiDbModeType _multiDbMode;
    int _numClosestHitsWanted;

    bool handle_d();
    bool handle_D();
    bool handle_io();
    bool handle_iu();
    bool handle_id();
    bool handle_fu();
    bool handle_fd();
    bool handle_N();
    bool handle_t();
    bool handle_mdb();
    bool handle_k();
};


#endif /* CONTEXTCLOSEST_H_ */
