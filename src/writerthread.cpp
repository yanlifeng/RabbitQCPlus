#include "writerthread.h"
#include "util.h"
#include <memory.h>
#include <unistd.h>
#include "read.h"

WriterThread::WriterThread(Options *opt, string filename) {
    mOptions = opt;

    mWriter1 = NULL;

    mInputCounter = 0;
    mOutputCounter = 0;
    mInputCompleted = false;
    mFilename = filename;

    mRingBuffer = new char *[PACK_NUM_LIMIT];
    memset(mRingBuffer, 0, sizeof(char *) * PACK_NUM_LIMIT);
    mRingBufferSizes = new size_t[PACK_NUM_LIMIT];
    memset(mRingBufferSizes, 0, sizeof(size_t) * PACK_NUM_LIMIT);
    initWriter(filename);
}

WriterThread::~WriterThread() {
    cleanup();
    delete mRingBuffer;
}

bool WriterThread::isCompleted() {
    return mInputCompleted && (mOutputCounter == mInputCounter);
}

bool WriterThread::setInputCompleted() {
    mInputCompleted = true;
    return true;
}

void WriterThread::output() {
    if (mOutputCounter >= mInputCounter) {
        usleep(100);
    }
    unsigned long nowTotle = 0;
    while (mOutputCounter < mInputCounter) {
        nowTotle += mRingBufferSizes[mOutputCounter];
        mWriter1->write(mRingBuffer[mOutputCounter], mRingBufferSizes[mOutputCounter]);
        delete[] mRingBuffer[mOutputCounter];
        mRingBuffer[mOutputCounter] = NULL;
        mOutputCounter++;
    }
//    printf("nowTotle %ld b\n", nowTotle);
}

void WriterThread::input(char *data, size_t size) {
    mRingBuffer[mInputCounter] = data;
    mRingBufferSizes[mInputCounter] = size;
    mInputCounter++;
}

//TODO
void WriterThread::input(vector<Read *> newOut) {
//
//    for (int i = 0; i < newOut.size(); i++) {
//        Read *now = newOut[i];
//        memcpy(nowPos, now->mName.c_str(), now->mName.size());
//        nowPos += now->mName.size();
//        *nowPos = '\n';
//        nowPos++;
//        memcpy(nowPos, now->mSeq.mStr.c_str(), now->mSeq.length());
//        nowPos += now->mSeq.length();
//        *nowPos = '\n';
//        nowPos++;
//        memcpy(nowPos, now->mStrand.c_str(), now->mStrand.size());
//        nowPos += now->mStrand.size();
//        *nowPos = '\n';
//        nowPos++;
//        memcpy(nowPos, now->mQuality.c_str(), now->mQuality.size());
//        nowPos += now->mQuality.size();
//        *nowPos = '\n';
//        nowPos++;
//        delete now;
//    }
//    mRingBuffer[mInputCounter] = data;
//    mRingBufferSizes[mInputCounter] = size;
//    mInputCounter++;
}

void WriterThread::cleanup() {
    deleteWriter();
}

void WriterThread::deleteWriter() {
    if (mWriter1 != NULL) {
        delete mWriter1;
        mWriter1 = NULL;
    }
}

void WriterThread::initWriter(string filename1) {
    deleteWriter();
    mWriter1 = new Writer(filename1, mOptions->compression);
}

void WriterThread::initWriter(ofstream *stream) {
    deleteWriter();
    mWriter1 = new Writer(stream);
}

void WriterThread::initWriter(gzFile gzfile) {
    deleteWriter();
    mWriter1 = new Writer(gzfile);
}

long WriterThread::bufferLength() {
    return mInputCounter - mOutputCounter;
}
