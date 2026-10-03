/*
 * PushBackStream.cpp
 *
 *  Created on: Mar 18, 2013
 *      Author: nek3d
 */
#include <cstdio>

#include "PushBackStreamBuf.h"
#include <cstring>
#include <iostream>

PushBackStreamBuf::PushBackStreamBuf(streambuf* primary_stream)
: streambuf(),
  _primary_stream(primary_stream)
 {

}

PushBackStreamBuf::~PushBackStreamBuf()
{
}

// Consume and return the next character, drawing from the pushed-back buffer
// first and falling through to the underlying stream once it is empty.
int PushBackStreamBuf::sbumpc()
 {
	int c = 0;
	if (!_buffer.empty()) {
		c = _buffer.front()->value();
		_buffer.pop_front();
		return c;
	}
	c = _primary_stream->sbumpc();
	return c;
}

// Replace the pushed-back buffer with newBuf; its contents are returned by
// subsequent reads before the underlying stream is touched again.
void PushBackStreamBuf::pushBack(const BTlist<int> &newBuf)
{
	_buffer = newBuf;
}

// Peek at the next character (without consuming it), preferring the pushed-back
// buffer and otherwise the underlying stream.
int PushBackStreamBuf::underflow()
{
	int c = 0;

	if(!_buffer.empty()) {
		c=_buffer.front()->value();
		return c;
	}
	c= _primary_stream->sgetc();
	return c;
}
