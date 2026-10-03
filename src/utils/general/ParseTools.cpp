#include "ParseTools.h"
#include <climits>
#include <cctype>
#include <cstring>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <sstream>


// Radix used by the hand-rolled integer parser in str2chrPos().
static const int DECIMAL_BASE = 10;


//This functions recognizes only numbers with digits, plus sign, minus sign, decimal point, e, or E. Hexadecimal and pointers not currently supported.
bool isNumeric(const string &str) {
	bool hasDigits = false;
	for (int i=0; i < (int)str.size(); i++) {
		char currChar = str[i];
		if (!(isdigit(currChar) || currChar == '-' || currChar == '.' || currChar == '+' || currChar == 'e' || currChar == 'E')) {
			return false;
		}
		hasDigits |= isdigit(currChar);
	}
	return hasDigits;
}

//As above, but does not allow decimal points (integers only, optionally signed
//or in scientific notation; at least one digit is required).
bool isInteger(const string &str) {
	bool hasDigits = false;
	for (int i=0; i < (int)str.size(); i++) {
		char currChar = str[i];
		if (!(isdigit(currChar) || currChar == '-' || currChar == '+' || currChar == 'e' || currChar == 'E')) {
			return false;
		}
		hasDigits |= isdigit(currChar);
	}
	return hasDigits;
}

// Strip any trailing carriage-return / line-feed characters from str in place.
void trimNewlines(string& str) {
	size_t len = str.length();
	while (len > 0 && (str[len-1] == '\n' || str[len-1] == '\r'))
		len--;
	str.resize(len);
}


// Convenience overload: parse an entire std::string as a chromosome position.
CHRPOS str2chrPos(const string &str) {
	return str2chrPos(str.c_str(), str.size());
}

// Fast base-10 parser for chromosome positions: reads an optional leading '-'
// then decimal digits, with a fallback to strtod for scientific notation.
// Exits with an error on any other trailing characters.
CHRPOS str2chrPos(const char * __restrict str, size_t ulen) {

	if (ulen == 0) {
		ulen = strlen(str);
	}

	const char* endpos = str;
	long long result = 0;
	bool neg = false;
	char last = 0;

	if(*endpos == '-') neg = true, endpos ++;

	for(;(last = *endpos); endpos ++) {
		if(last < '0' || last > '9') break;
		result = result * DECIMAL_BASE + last - '0';
	}

	if(last) {
		if(*endpos == 'e' || *endpos == 'E') {
			char* endpos = NULL;
			CHRPOS ret = (CHRPOS)strtod(str, &endpos);

			if(endpos && *endpos == 0) {
				return ret;
			}
		}
		fprintf(stderr, "***** ERROR: illegal number \"%s\". Exiting...\n", str);
		exit(1);
	}

	return neg?-result:result;
}

// Build a string whose characters are the (truncated) integer values in vec.
string vectorIntToStr(const vector<int> &vec) {
	string str;
	str.reserve(vec.size());
	for (int i=0; i < (int)vec.size(); i++) {
		str += (char)(vec[i]);
	}
	return str;
}


// Return true if the line looks like a header/comment rather than a data row:
// it starts with '>', '!' or '#', or begins with one of the UCSC header
// keywords ("chrom"/"chr" followed by whitespace, "browser", "track",
// "visibility"). Two implementations exist; the x86 variant hand-unrolls the
// keyword comparison for speed but is functionally identical to the portable one.
#if defined(__i386__) || defined(__x86_64__)
bool isHeaderLine(const string &line) {
	if (line[0] == '>') {
		return true;
	}
	if (line[0] == '!') {
		return true;
	}
	if (line[0] == '#') {
		return true;
	}
	
	if(line.length() > 4) {
		uint32_t peek = *(uint32_t*)line.c_str() | 0x20202020u;
		const char* full_text = NULL;
		bool require_space = false;
		bool require_digit = false;
		switch(peek) {
			case 0x6f726863: 
				full_text = "chrom";
				require_space = true;
				require_digit = false;
				break;
			case 0x20726863:
			case 0x09726863:
				full_text = "chr";
				require_space = require_digit = true;
				break;
			case 0x776f7262:
				full_text = "browser";
				break;
			case 0x63617274:
				full_text = "track";
				break;
			case 0x69736976:
				full_text = "visibility";
				break;
			default:
				return false;
		}
		if(full_text) {
			const char* ptr = NULL;
			for(ptr = line.c_str(); *ptr && *full_text; ptr ++, full_text ++) {
				char c = *ptr;
				if(c >= 'A' && c <= 'Z') c += 32;
				if(c != *full_text) return false;
			}
			if(require_space && !isspace(*(ptr++))) return false;
			if(require_digit && !isdigit(*(ptr++))) return false;
			return true;
		}
	}

	return false;
}
#else
bool isHeaderLine(const string &line) {
	if (line[0] == '>') {
		return true;
	}
	if (line[0] == '!') {
		return true;
	}
	if (line[0] == '#') {
		return true;
	}
	
	string tmp = line;
	transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
	//allow chr chrom to start a header line
	if (memcmp(tmp.c_str(), "chrom", 5) == 0 && isspace(tmp[5]) && ! isdigit(tmp[6])) {
		return true;
	}
	//allow chr chrom to start a header line
	if (memcmp(tmp.c_str(), "chr", 3) == 0 && isspace(tmp[3]) && ! isdigit(tmp[4])) {
		return true;
	}
	//UCSC file headers can also start with the words "browser" or "track", followed by a whitespace character.
	if (memcmp(tmp.c_str(), "browser", 7) == 0) {
		return true;
	}
	if (memcmp(tmp.c_str(), "track", 5) == 0) {
		return true;
	}
	if (memcmp(tmp.c_str(), "visibility", 10) == 0) {
		return true;
	}
	return false;
}
#endif
