// MIT License
//
// Copyright (c) 2017
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//
//
// Created by kellerberrin on 2/09/17.
//

#ifndef READSAMFILE_THREAD_H
#define READSAMFILE_THREAD_H

#endif //READSAMFILE_THREAD_H


#include <string>
#include <vector>
#include <queue>


class QueueMessage {

public:

    QueueMessage() = default;
    virtual ~QueueMessage() = default;

    virtual bool isEOF() = 0;

};




class EOFMessage : public QueueMessage {

public:

    EOFMessage() = default;
    virtual ~EOFMessage() = default;

    inline bool isEOF() { return true; }

};




class SamRecord : public QueueMessage {

public:

    SamRecord(const std::string& samRecord) { _parseSamRecord(samRecord); }
    virtual ~SamRecord() = default;

    const std::string qName;
    const std::string rName;
    const std::string pos;
    const std::string mapQuality;
    const std::string cigar;
    const std::string rNext;
    const std::string pNext;
    const std::string tLen;
    const std::string sequence;
    const std::string quality;
    const std::vector<const std::string> optFlags;

    inline bool isEOF() { return false; }

private:

    void _parseSamRecord(const std::string& samRecord) throw;

};







class IOProducerConsumerQueue {

public:

    IOProducerConsumerQueue(unsigned long HighWaterMark, unsigned long LowWaterMark):
            _highWaterMark(HighWaterMark), _lowWaterMark(LowWaterMark) {}
    virtual ~IOProducerConsumerQueue() = default;

    void enQueueMessage(const QueueMessage& message);     // Atomic. Blocks if queue size == high water mark
    const QueueMessage& deQueueMessage();      // Atomic. Blocks if the queue is empty.

private:

    unsigned long _highWaterMark;
    unsigned long _lowWaterMark;
    std::queue<const QueueMessage&> _messageQueue;

};


