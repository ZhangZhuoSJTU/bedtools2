/*****************************************************************************
  VecOps.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef VECTOROPS_H
#define VECTOROPS_H

#include <vector>
#include <map>
#include <math.h>
#include <stdint.h>
#include <string>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <limits.h>

using namespace std;

typedef int64_t CHRPOS;

//************************************************
// VectorOps Class methods and elements
//************************************************
//
// VectorOps wraps a vector of string tokens (typically one column extracted
// from a set of records) and exposes a family of summary operations over it:
// numeric reductions (sum, mean, stddev, median, min/max), categorical
// reductions (mode, anti-mode, distinct, counts) and textual joins (collapse,
// concat, frequency histograms).
//
// Numeric operations lazily materialise a parallel vector<double> (_vecd) the
// first time one is requested; several of the methods mutate the underlying
// string vector (e.g. by sorting it), so an instance is intended to be used
// for a single pass of reductions rather than reused across unrelated calls.
//************************************************
class VectorOps {

public:

    // Constructor. Initialize with a vector of strings.
    VectorOps(const vector<string> &vec);

    // Destructor
    ~VectorOps(void);

    // user-interface
    void setDelim(const string &str) { _delimStr = str; }
    
    // return the total of the values in the vector
    double GetSum(void);             
    // return the average value in the vector
    double GetMean(void);             
     // return the standard deviation
    double GetStddev(void);          
    // return the sample standard deviation 
    double GetSstddev(void);          
    // return the median value in the vector
    double GetMedian(void);           
    // return the most common value in the vector
    string GetMode(void);             
    // return the least common value in the vector
    string GetAntiMode(void);         
    // return the minimum element of the vector
    double GetMin(void);              
    // return the maximum element of the vector
    double GetMax(void);              
    // return the minimum absolute value of the vector
    double GetAbsMin(void);              
    // return the maximum absolute value of the vector
    double GetAbsMax(void);              
    // return the count of element in the vector
    uint32_t GetCount(void);          
    // return a the count of _unique_ elements in the vector
    uint32_t GetCountDistinct(void);  
    // return a delimiter-separated list of elements
    string GetCollapse(string delimiter = ",");
    // return a concatenation of all elements in the vector
    string GetConcat(void);           
    // return a comma-separated list of the _unique_ elements
    string GetDistinct(void);
    // return a comma-separated list of the _unique_ elements, sorted numerically
    string GetDistinctSortNum(bool asc = true);
    // return a histogram of values and their freqs. in desc. order of frequency
    string GetFreqDesc(void);         
    // return a histogram of values and their freqs. in asc. order of frequency
    string GetFreqAsc(void);          
    // return the first value in the list
    string GetFirst(void);
    // return the last value in the list
    string GetLast(void);
    
private:
    vector<string> _vecs;   // the raw string tokens being summarised
    vector<double> _vecd;    // cached numeric view of _vecs, filled on demand
    uint32_t _size;          // number of tokens supplied at construction
    string _delimStr;        // delimiter used when joining elements into output
};

#endif /* VECTOROPS_H */
