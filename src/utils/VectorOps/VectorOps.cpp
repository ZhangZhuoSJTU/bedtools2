/*****************************************************************************
  VectorOps.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Purpose: Container to perform various useful operations 
           (e.g., min, max, count, distinct, etc.) on vectors.

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#include "VectorOps.h"

// Process-exit status used when a fatal, unrecoverable error occurs (e.g. a
// value in the vector cannot be parsed as a number). Replaces the bare
// literal "1" that was previously scattered through the conversion functors.
static const int VECTOROPS_ERROR_EXIT_CODE = 1;

// Character written between a value and its observed frequency when emitting
// a histogram (e.g. "valueA:3,valueB:1"). Previously a bare ":" literal.
static const char FREQUENCY_SEPARATOR = ':';

// functor to convert a single string to a double.
// Parses `element` as a floating-point value; aborts the program with a
// diagnostic if the text is not a valid number.
double MakeDouble(const string &element) {
    std::istringstream valueStream(element);
    double parsedValue;
    if (!(valueStream >> parsedValue)) {
        cerr << "Error: Could not properly convert string to numeric (\"" + element + "\")" << endl;
        exit(VECTOROPS_ERROR_EXIT_CODE);
    }
    return parsedValue;
}

// functor to convert the abs() of a single string to a double.
// Behaves like MakeDouble, but returns the absolute value of the parsed
// number. Used by the "absolute" min/max reductions.
double MakeAbsDouble(const string &element) {
    std::istringstream valueStream(element);
    double parsedValue;
    if (!(valueStream >> parsedValue)) {
        cerr << "Error: Could not properly convert string to numeric (\"" + element + "\")" << endl;
        exit(VECTOROPS_ERROR_EXIT_CODE);
    }
    return fabs(parsedValue);
}

// Comparator that orders (frequency, value) pairs by descending frequency.
// Used to sort a histogram so the most common values appear first.
struct ValueGreaterThan
{
    bool operator()( const vector< pair<int, string> >::value_type& lhs,
        const vector< pair<int, string> >::value_type& rhs ) const
    {
        return lhs.first > rhs.first;
    }
};

// Comparator that orders (frequency, value) pairs by ascending frequency.
// Used to sort a histogram so the least common values appear first.
struct ValueLessThan
{
    bool operator()( const vector< pair<int, string> >::value_type& lhs,
        const vector< pair<int, string> >::value_type& rhs ) const
    {
        return lhs.first < rhs.first;
    }
};

// Constructor
VectorOps::VectorOps(const vector<string> &vec)
: _vecs(vec)
, _vecd()
, _size(vec.size())
, _delimStr(",")
{
    _vecd.reserve(vec.size());
}

// Destructor
VectorOps::~VectorOps(void) {
}

double VectorOps::GetSum(void) 
{
    // convert the vec of strings to a vec of doubles
    transform(_vecs.begin(), _vecs.end(), back_inserter(_vecd), MakeDouble);
    return accumulate(_vecd.begin(), _vecd.end(), 0.0);
}

double VectorOps::GetMean(void) 
{
    // convert the vec of strings to a vec of doubles
    transform(_vecs.begin(), _vecs.end(), back_inserter(_vecd), MakeDouble);
    return accumulate(_vecd.begin(), _vecd.end(), 0.0) / _size;
}

// Population standard deviation: sqrt of the mean squared deviation from the
// mean, dividing the summed squared deviations by N.
double VectorOps::GetStddev(void)
{
    const double mean = GetMean();
    // Accumulate the sum of squared deviations from the mean.
    double sumSquaredDeviations = 0.0;
    vector<double>::const_iterator valueIter = _vecd.begin();
    vector<double>::const_iterator valuesEnd = _vecd.end();
    for (; valueIter != valuesEnd; ++valueIter) {
        sumSquaredDeviations += pow((*valueIter - mean), 2);
    }
    const double variance = sumSquaredDeviations / _vecd.size();
    return sqrt(variance);
}

// Sample standard deviation: identical to GetStddev, except the summed
// squared deviations are divided by (N - 1) (Bessel's correction).
double VectorOps::GetSstddev(void)
{
    const double mean = GetMean();
    // Accumulate the sum of squared deviations from the mean.
    double sumSquaredDeviations = 0.0;
    vector<double>::const_iterator valueIter = _vecd.begin();
    vector<double>::const_iterator valuesEnd = _vecd.end();
    for (; valueIter != valuesEnd; ++valueIter) {
        sumSquaredDeviations += pow((*valueIter - mean), 2);
    }
    const double variance = sumSquaredDeviations / (_vecd.size() - 1);
    return sqrt(variance);
}

double VectorOps::GetMedian(void) 
{
    // convert the vec of strings to a vec of doubles
    transform(_vecs.begin(), _vecs.end(), back_inserter(_vecd), MakeDouble);
    
    double median = 0.0;
    sort(_vecd.begin(), _vecd.end());
    CHRPOS numValues = _vecd.size();
    if ((numValues % 2) > 0) {
        // Odd number of values: the median is the single middle element.
        CHRPOS mid;
        mid = numValues / 2;
        median = _vecd[mid];
    }
    else {
        // Even number of values: the median is the mean of the two middle
        // elements.
        CHRPOS midLow, midHigh;
        midLow = (numValues / 2) - 1;
        midHigh = (numValues / 2);
        median = (_vecd[midLow] + _vecd[midHigh]) / 2.0;
    }
    return median;
}

double VectorOps::GetMin(void)
{
    // convert the vec of strings to a vec of doubles
    transform(_vecs.begin(), _vecs.end(), back_inserter(_vecd), MakeDouble);
    return *min_element( _vecd.begin(), _vecd.end() );
}

double VectorOps::GetAbsMin(void)
{
    // convert the vec of strings to a vec of doubles
    transform(_vecs.begin(), _vecs.end(), back_inserter(_vecd), MakeAbsDouble);
    return *min_element( _vecd.begin(), _vecd.end() );
}

double VectorOps::GetMax(void)
{
    // convert the vec of strings to a vec of doubles
    transform(_vecs.begin(), _vecs.end(), back_inserter(_vecd), MakeDouble);
    return *max_element( _vecd.begin(), _vecd.end() );
}

double VectorOps::GetAbsMax(void)
{
    // convert the vec of strings to a vec of doubles
    transform(_vecs.begin(), _vecs.end(), back_inserter(_vecd), MakeAbsDouble);
    return *max_element( _vecd.begin(), _vecd.end() );
}

// Return the most frequently occurring value (the statistical mode).
// The values are sorted so equal entries are adjacent, then we walk the
// sorted list counting run lengths and remembering the longest run seen.
string VectorOps::GetMode(void)
{
    string mode;
    int maxRunLength = 0;
    int currRunLength = 0;

    sort(_vecs.begin(), _vecs.end());

    string prevValue = "";
    vector<string>::const_iterator valueIter = _vecs.begin();
    while (valueIter != _vecs.end())
    {
        if (*valueIter != prevValue && prevValue != "")
        {
            // We just finished a run of identical values; see if it is the
            // longest run encountered so far.
            if (currRunLength > maxRunLength)
            {
                maxRunLength = currRunLength;
                mode = prevValue;
            }
            currRunLength = 1;
        }
        else { currRunLength++; }
        prevValue = *valueIter;
        ++valueIter;
    }
    // Account for the final run, which ends when the loop exits.
    if (currRunLength > maxRunLength)
        mode = prevValue;
    return mode;
}


// Return the least frequently occurring value (the "anti-mode").
// Mirrors GetMode, but tracks the shortest run of identical values instead
// of the longest.
string VectorOps::GetAntiMode(void)
{
    string antimode;
    int minRunLength = INT_MAX;
    int currRunLength = 0;

    sort(_vecs.begin(), _vecs.end());

    string prevValue = "";
    vector<string>::const_iterator valueIter = _vecs.begin();
    while (valueIter != _vecs.end())
    {
        if (*valueIter != prevValue && prevValue != "")
        {
            // We just finished a run of identical values; see if it is the
            // shortest run encountered so far.
            if (currRunLength < minRunLength)
            {
                minRunLength = currRunLength;
                antimode = prevValue;
            }
            currRunLength = 1;
        }
        else { currRunLength++; }
        prevValue = *valueIter;
        ++valueIter;
    }
    // Account for the final run, which ends when the loop exits.
    if (currRunLength < minRunLength)
        antimode = prevValue;
    return antimode;
}

uint32_t VectorOps::GetCount(void)
{
    return _vecs.size();
}

uint32_t VectorOps::GetCountDistinct(void)
{
    sort( _vecs.begin(), _vecs.end() );
    _vecs.erase( unique( _vecs.begin(), _vecs.end() ), _vecs.end() );
    return _vecs.size();
}

string VectorOps::GetCollapse(string delimiter)
{
    ostringstream collapse;
    for( size_t i = 0; i < _vecs.size(); i++ ) {
        if (i>0)
            collapse << delimiter;
        collapse << _vecs[i];
    }
    return collapse.str();
}

string VectorOps::GetConcat(void)
{
    ostringstream collapse;
    for( size_t i = 0; i < _vecs.size(); i++ )
        collapse << _vecs[i];
    return collapse.str();
}

string VectorOps::GetDistinct(void)
{
    ostringstream distinct;
    // remove duplicate entries from the vector
    // http://stackoverflow.com/questions/1041620/most-efficient-way-to-erase-duplicates-and-sort-a-c-vector
    sort( _vecs.begin(), _vecs.end() );
    _vecs.erase( unique( _vecs.begin(), _vecs.end() ), _vecs.end() );
    
    for( size_t i = 0; i < _vecs.size(); i++ ) {
        if (i>0)
            distinct << _delimStr;
        distinct << _vecs[i];
    }
    return distinct.str();
}

string VectorOps::GetDistinctSortNum(bool ascending)
{
    ostringstream distinct;
    // remove duplicate entries from the vector
    // http://stackoverflow.com/questions/1041620/most-efficient-way-to-erase-duplicates-and-sort-a-c-vector
    // convert the vec of strings to a vec of doubles
    transform(_vecs.begin(), _vecs.end(), back_inserter(_vecd), MakeDouble);

    if (ascending) {
    	sort( _vecd.begin(), _vecd.end(), less<double>());
    } else {
    	sort( _vecd.begin(), _vecd.end(), greater<double>());

    }
    _vecd.erase( unique( _vecd.begin(), _vecd.end() ), _vecd.end() );

    for( size_t i = 0; i < _vecd.size(); i++ ) {
        if (i>0)
            distinct << _delimStr;
        distinct << _vecd[i];
    }
    return distinct.str();
}

string VectorOps::GetFreqDesc(void)
{
    // compute the frequency of each unique value
    map<string, int> freqs;
    vector<string>::const_iterator dIt  = _vecs.begin();
    vector<string>::const_iterator dEnd = _vecs.end();
    for (; dIt != dEnd; ++dIt) {
        freqs[*dIt]++;
    }
    // pair for the num times a values was
    // observed (1) and the value itself (2)
    pair<int, string> freqPair;
    vector< pair<int, string> > freqList;

    // create a list of pairs of all the observed values (second)
    // and their occurences (first)
    map<string,int>::const_iterator mapIter = freqs.begin();
    map<string,int>::const_iterator mapEnd  = freqs.end();
    for(; mapIter != mapEnd; ++mapIter)
        freqList.push_back( make_pair(mapIter->second, mapIter->first) );

    // sort the list of pairs in the requested order by the frequency
    // this will make the value that was observed least/most bubble to the top
    sort(freqList.begin(), freqList.end(), ValueGreaterThan());

    // record all of the values and their frequencies.
    ostringstream buffer;
    vector< pair<int, string> >::const_iterator iter    = freqList.begin();
    vector< pair<int, string> >::const_iterator iterEnd = freqList.end();
    for (; iter != iterEnd; ++iter)
        buffer << iter->second << FREQUENCY_SEPARATOR << iter->first << _delimStr;

    return buffer.str();
}


string VectorOps::GetFreqAsc(void)
{
    // compute the frequency of each unique value
    map<string, int> freqs;
    vector<string>::const_iterator dIt  = _vecs.begin();
    vector<string>::const_iterator dEnd = _vecs.end();
    for (; dIt != dEnd; ++dIt) {
        freqs[*dIt]++;
    }
    // pair for the num times a values was
    // observed (1) and the value itself (2)
    pair<int, string> freqPair;
    vector< pair<int, string> > freqList;

    // create a list of pairs of all the observed values (second)
    // and their occurences (first)
    map<string,int>::const_iterator mapIter = freqs.begin();
    map<string,int>::const_iterator mapEnd  = freqs.end();
    for(; mapIter != mapEnd; ++mapIter)
        freqList.push_back( make_pair(mapIter->second, mapIter->first) );

    // sort the list of pairs in the requested order by the frequency
    // this will make the value that was observed least/most bubble to the top
    sort(freqList.begin(), freqList.end(), ValueLessThan());

    // record all of the values and their frequencies.
    ostringstream buffer;
    vector< pair<int, string> >::const_iterator iter    = freqList.begin();
    vector< pair<int, string> >::const_iterator iterEnd = freqList.end();
    for (; iter != iterEnd; ++iter)
        buffer << iter->second << FREQUENCY_SEPARATOR << iter->first << _delimStr;

    return buffer.str();
}


string VectorOps::GetFirst(void)
{
    return _vecs[0];
}


string VectorOps::GetLast(void)
{
    return _vecs[_vecs.size() - 1];
}
