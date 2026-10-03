/*
 * Tokenizer.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: nek3d
 */
#include "Tokenizer.h"
#include <cstring>
#include <cstdio>
#include <iostream>
Tokenizer::Tokenizer()
: _numExpectedElems(0),
  _keepFinalIncElem(USE_NOW) {
}

Tokenizer::~Tokenizer() {
}

// Hint the expected number of fields so buffers can be sized appropriately.
void Tokenizer::setNumExpectedItems(int newSize) {
	_numExpectedElems = newSize;
}

// Count how many fields str would split into on the given delimiter, without
// disturbing the tokenizer's stored elements.
int Tokenizer::getNumFields(const string &str, char delimiter)
{
    // http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c/236803#236803
    // NOTE: this approach intentionally allows consecutive delimiters
    std::stringstream ss(str);
    std::string item;
    vector<string> tmp;
    while(getline(ss, item, delimiter)) {
        tmp.push_back(item);  
    }
    return tmp.size();
}

// Split str on the delimiter into the stored element list and return the field
// count. For compressed input split on newlines, a trailing partial line (one
// not terminated by '\n') is dropped so only complete lines are sampled.
int Tokenizer::tokenize(const string &str, char delimiter, bool eofHit, bool isCompressed) {

    // http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c/236803#236803
    // NOTE: this approach intentionally allows consecutive delimiters
    std::stringstream ss(str);
    std::string item;
    _elems.clear();
    while(getline(ss, item, delimiter)) {
    	// when dealing with compressed input, we can sometimes token lines such
    	// that we capture a partial line (i.e., it didn't end in '\n'). Only
    	// sample the complete lines.
    	if (isCompressed)
    	{
    		if (delimiter == '\n' && ss.eof() == false)
    		{  
    			_elems.push_back(item);
    		}
    		else if (delimiter != '\n')
    		{
    			_elems.push_back(item);
    		}
    	}
    	else
    	{
    		_elems.push_back(item);
    	}
    }
    return _elems.size();


	// int strLen = (int)str.size();

	// int startPos = 0;
	// int currPos = 0;
	// int currIdx = 0;
	// while (startPos < strLen) {
	// 	while (str[currPos] != delimiter && currPos < strLen) {
	// 		currPos++;
	// 	}
	// 	if (currPos > startPos) {
	// 		if ((currPos == strLen && _keepFinalIncElem != USE_NOW) &&
	// 				( (!(delimiter == '\n' && eofHit)) || isCompressed)) {
	// 			//we found an incomplete final element.
	// 			// if we're ignoring incomplete elems, do nothing with it.
	// 			currIdx--; //make sure it's not included in the final count of valid elems.

	// 		} else {
	// 			string *newStr = fetchElem(currIdx);
	// 			newStr->assign(str.c_str() + startPos, min(currPos, strLen) - startPos);

	// 			// If splitting lines, strip any white space from the end of the line
	// 			// including DOS newline characters and excess tabs.
	// 			if (delimiter == '\n') {
	// 				int lastPos = newStr->size();
	// 				while (isspace(newStr->at(lastPos-1))) 
	// 				{
	// 					lastPos--;
	// 				}
	// 				newStr->resize(lastPos);
	// 			}
	// 		}
	// 	}
	// 	startPos = currPos +1;
	// 	currPos = startPos;
	// 	currIdx++;
	// }
	// _numValidElems = currIdx;
	// return currIdx;
}

// Choose whether a trailing element that ended without a delimiter is kept.
void Tokenizer::setKeepFinalIncompleteElem(lastElemCode code) {
	_keepFinalIncElem = code;
}

// Return the element at the given index (no bounds checking).
string Tokenizer::fetchElem(int idx)
{
	return _elems[idx];
}

