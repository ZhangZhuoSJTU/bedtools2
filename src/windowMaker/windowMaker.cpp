/*****************************************************************************
windowMaker.cpp

(c) 2009 - Aaron Quinlan
Hall Laboratory
Department of Biochemistry and Molecular Genetics
University of Virginia
aaronquinlan@gmail.com

Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "windowMaker.h"

using namespace WindowMakerConstants;

// Fixed-size constructor: windows of width "size" advanced by "step".
WindowMaker::WindowMaker(string &fileName, ID_METHOD id_method, INPUT_FILE_TYPE input_file_type, uint32_t size, uint32_t step, bool reverse)
: _size(size)
, _step(step)
, _count(0)
, _reverse(reverse)
, _window_method(FIXED_WINDOW_SIZE)
, _id_method(id_method)
{
    if (input_file_type==GENOME_FILE)
        MakeWindowsFromGenome(fileName);
    else
        MakeWindowsFromBED(fileName);
}

// Fixed-count constructor: split each interval into "count" windows.
WindowMaker::WindowMaker(string &fileName, ID_METHOD id_method, INPUT_FILE_TYPE input_file_type, uint32_t count, bool reverse)
    : _size(0)
    , _step(0)
    , _count(count)
    , _reverse(reverse)
    , _window_method(FIXED_WINDOW_COUNT)
    , _id_method(id_method)
{
    if (input_file_type==GENOME_FILE)
        MakeWindowsFromGenome(fileName);
    else
        MakeWindowsFromBED(fileName);
}

WindowMaker::~WindowMaker(void) {}


// Treat every chromosome of the genome file as a single interval
// spanning [0, chromSize) and window each one.
void WindowMaker::MakeWindowsFromGenome(const string& genomeFileName) {

    NewGenomeFile *_genome = new NewGenomeFile(genomeFileName);

    // get a list of the chroms in the user's genome
    vector<string> chromList = _genome->getChromList();

    // process each chrom in the genome
    for (size_t c = 0; c < chromList.size(); ++c) {
        string chrom = chromList[c];

        BED bed(chrom,0,_genome->getChromSize(chrom),chrom,"","");
        MakeBEDWindow(bed);
    }
}

// Window each valid feature read from the input BED file.
void WindowMaker::MakeWindowsFromBED(string& bedFileName) {
    BedFile bf(bedFileName);
    bf.Open();

    BED bed;
    while (bf.GetNextBed(bed)) {
        if (bf._status == BED_VALID)
            MakeBEDWindow(bed);
    }
    bf.Close();
}

// Dispatch to the fixed-size or fixed-count strategy chosen at construction.
void WindowMaker::MakeBEDWindow(const BED& interval)
{
    if (_window_method==FIXED_WINDOW_SIZE)
        MakeFixedSizeWindow(interval);
    else
        MakeFixedCountWindow(interval);
}

// Compute how many fixed-size windows (given the current _step) are needed
// to cover the interval, rounding up so a short trailing window still counts.
uint32_t WindowMaker::CalculateWindows(const BED& interval) {
    uint32_t num_windows = (interval.end - interval.start) / _step;
    if ((interval.end - interval.start) % _step > 0) {
        // add 1 to num_windows if the last window is less than _step
        num_windows += 1;
    }
    return num_windows;
}
// Emit fixed-width windows, stepping across the interval by _step.
void WindowMaker::MakeFixedSizeWindow(const BED& interval) {
    uint32_t i = FIRST_WINDOW_INDEX;
    uint32_t num_windows = CalculateWindows(interval);
    for (uint32_t start = interval.start; start <= interval.end; start += _step, ++i) {
        string name = GenerateID(interval, i, num_windows, _reverse);
        if ((start + _size) <= interval.end) {
            cout << interval.chrom << "\t" << start << "\t" << start + _size << name << endl;
        }
        else if (start < interval.end) {
            cout << interval.chrom << "\t" << start << "\t" << interval.end << name << endl;
        }
    }
}

// Partition the interval into exactly _count equal windows, extending the
// final window to absorb any remainder so the whole interval is covered.
void WindowMaker::MakeFixedCountWindow(const BED& interval) {
    CHRPOS interval_size = interval.end - interval.start ;
    CHRPOS window_size = (interval_size)/_count; // integer version of ceil(interval_size/_count)

    if (window_size == 0)
    {
        cerr << "WARNING: Interval " 
             << interval.chrom << ":" 
             << interval.start << "-" 
             << interval.end
             << " is smaller than the number of windows requested. Skipping."
             << endl;
        return;
    } 
    if (interval_size == 0)
    {
        cerr << "WARNING: Interval " 
             << interval.chrom << ":" 
             << interval.start << "-" 
             << interval.end
             << " has zero length and cannot be partitioned. Skipping."
             << endl;
        return;
    }


    uint32_t i = FIRST_WINDOW_INDEX;
    for (CHRPOS start = interval.start; start < interval.end; start += window_size, ++i) {
        string name = GenerateID(interval, i, _count, _reverse);
        CHRPOS end = min(start + window_size,interval.end);

        // extend range of last interval if necessary
        if (i == _count) 
        {
            end = interval.end;
            cout << interval.chrom << "\t" << start << "\t" << end << name << endl;
            break;
        }
        cout << interval.chrom << "\t" << start << "\t" << end << name << endl;
    }
}

// Build the optional trailing ID column for a window according to the
// user's chosen ID_METHOD (none, window number, source id, or both).
// When _reverse is requested, window numbers count down instead of up.
string WindowMaker::GenerateID(const BED& interval, uint32_t window_index, uint32_t num_windows, bool _reverse) const {
    stringstream s;
    switch(_id_method) {
    case ID_SOURCE_ID:
         s << "\t" << interval.name;
         break;
    case ID_WINDOW_NUMBER:
         if (_reverse == true && num_windows > 0) {
            s << "\t" << num_windows - window_index + 1;
         } else {
            s << "\t" << window_index;
         }
         break;
    case ID_SOURCE_ID_WINDOW_NUMBER:
         if (_reverse == true && num_windows > 0) {
            s << "\t" << interval.name << "_" << num_windows - window_index + 1;
         } else {
            s << "\t" << interval.name << "_" << window_index;
         }
    default:
    case ID_NONE:
         break;
    }
    return s.str();
}
