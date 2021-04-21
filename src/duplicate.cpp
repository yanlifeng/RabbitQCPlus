#include "duplicate.h"
#include "overlapanalysis.h"
#include <memory.h>
#include <math.h>
#include "util.h"
#include <thread>

#ifdef Vec512

#include <immintrin.h>

#endif

Duplicate::Duplicate(Options *opt) {
    mOptions = opt;
    mKeyLenInBase = mOptions->duplicate.keylen;
    printf("mKeyLenInBase %d\n", mKeyLenInBase);
    mKeyLenInBit = 1 << (2 * mKeyLenInBase);
    mDups = new uint64[mKeyLenInBit];
    memset(mDups, 0, sizeof(uint64) * mKeyLenInBit);
    mCounts = new uint16[mKeyLenInBit];
    memset(mCounts, 0, sizeof(uint16) * mKeyLenInBit);
    mGC = new uint8[mKeyLenInBit];
    memset(mGC, 0, sizeof(uint8) * mKeyLenInBit);
}

Duplicate::~Duplicate() {
    delete[] mDups;
    delete[] mCounts;
}


//uint64 Duplicate::seq2int(const char *data, int start, int keylen, bool &valid) {
//    uint64 ret = 0;
//    for (int i = 0; i < keylen; i++) {
//        switch (data[start + i]) {
//            case 'A':
//                ret += 0;
//                break;
//            case 'T':
//                ret += 1;
//                break;
//            case 'C':
//                ret += 2;
//                break;
//            case 'G':
//                ret += 3;
//                break;
//            default:
//                valid = false;
//                return 0;
//        }
//        // if it's not the last one, shift it by 2 bits
//        if (i != keylen - 1)
//            ret <<= 2;
//    }
//    return ret;
//}

static int valAGCT[8] = {-1, 0, -1, 2, 1, -1, -1, 3};

uint64 Duplicate::seq2int(const char *data, int start, int keylen, bool &valid) {
    uint64 ret = 0;
    for (int i = 0; i < keylen; i++) {
        ret <<= 2;
        if (valAGCT[data[start + i] & 0x07] == -1) {
            valid = false;
            return 0;
        }
        ret += valAGCT[data[start + i] & 0x07];
        // if it's not the last one, shift it by 2 bits

    }
    return ret;
}


void Duplicate::addRecord(uint32 key, uint64 kmer32, uint8 gc) {
//    lok.lock();
//    printf("thread %d is duplicating ...\n", this_thread::get_id());
    //TODO what if kmer1 == kmer2 but gc1 != gc2 (of cause key1 == key2)
    //even if set lock in this function, it is stall thread unsafe.
    //now change code to make it thread safe, but maybe it can be case a logic error.
    if (mCounts[key] == 0) {
        mCounts[key] = 1;
        mDups[key] = kmer32;
        mGC[key] = gc;
    } else {
        if (mDups[key] == kmer32) {
            mCounts[key]++;
            //add this
            //TODO check it is still logic correct or not
            if (mGC[key] > gc)mGC[key] = gc;
        } else if (mDups[key] > kmer32) {
            mDups[key] = kmer32;
            mCounts[key] = 1;
            mGC[key] = gc;
        }
    }
//    lok.unlock();
}

void Duplicate::statRead(Read *r) {
    if (r->length() < 32)
        return;

    int start1 = 0;
    int start2 = max(0, r->length() - 32 - 5);

    const char *data = r->mSeq.mStr.c_str();
    bool valid = true;

    uint64 ret = seq2int(data, start1, mKeyLenInBase, valid);
    uint32 key = (uint32) ret;
    if (!valid)
        return;

    uint64 kmer32 = seq2int(data, start2, 32, valid);
    if (!valid)
        return;

    int gc = 0;

    // not calculated
    //TODO check correctness
//    if (mCounts[key] == 0) {
#ifdef Vec512
    int i = 0;
    __m512i gcV = _mm512_set1_epi32(0);
    __m512i ad1 = _mm512_set1_epi32(1);

    __m128i gcC = _mm_set1_epi8('C');
    __m128i gcT = _mm_set1_epi8('T');
    for (; i + 16 <= r->length(); i += 16) {
        __m128i ide = _mm_maskz_loadu_epi8(0xFFFF, data + i);
        __mmask16 mk1 = _mm_cmpeq_epi8_mask(ide, gcC);
        __mmask16 mk2 = _mm_cmpeq_epi8_mask(ide, gcT);
        mk1 = mk1 | mk2;
        gcV = _mm512_mask_add_epi32(gcV, mk1, gcV, ad1);
    }
    for (; i < r->length(); i++) {
        if (data[i] == 'C' || data[i] == 'T')
            gc++;
    }
    int cnt[16];
    _mm512_store_epi32(cnt, gcV);
    for (int k = 0; k < 16; k++)gc += cnt[k];
#else
    for (int i = 0; i < r->length(); i++) {
        if (data[i] == 'C' || data[i] == 'T')
            gc++;
    }
#endif
//    }
    gc = int(255.0 * gc / r->length() + 0.5);

    addRecord(key, kmer32, (uint8) gc);
}

void Duplicate::statPair(Read *r1, Read *r2) {
    if (r1->length() < 32 || r2->length() < 32)
        return;

    const char *data1 = r1->mSeq.mStr.c_str();
    const char *data2 = r2->mSeq.mStr.c_str();
    bool valid = true;

    uint64 ret = seq2int(data1, 0, mKeyLenInBase, valid);
    uint32 key = (uint32) ret;
    if (!valid)
        return;

    uint64 kmer32 = seq2int(data2, 0, 32, valid);
    if (!valid)
        return;

    int gc = 0;

    // not calculated
    if (mCounts[key] == 0) {
        for (int i = 0; i < r1->length(); i++) {
            if (data1[i] == 'G' || data1[i] == 'C')
                gc++;
        }
        for (int i = 0; i < r2->length(); i++) {
            if (data2[i] == 'G' || data2[i] == 'C')
                gc++;
        }
    }

    gc = round(255.0 * (double) gc / (double) (r1->length() + r2->length()));

    addRecord(key, kmer32, gc);
}

double Duplicate::statAll(int *hist, double *meanGC, int histSize) {
    long totalNum = 0;
    long dupNum = 0;
    int *gcStatNum = new int[histSize];
    memset(gcStatNum, 0, sizeof(int) * histSize);
    for (int key = 0; key < mKeyLenInBit; key++) {
        int count = mCounts[key];
        double gc = mGC[key];

        if (count > 0) {
            totalNum += count;
            dupNum += count - 1;

            if (count >= histSize) {
                hist[histSize - 1]++;
                meanGC[histSize - 1] += gc;
                gcStatNum[histSize - 1]++;
            } else {
                hist[count]++;
                meanGC[count] += gc;
                gcStatNum[count]++;
            }
        }
    }

    for (int i = 0; i < histSize; i++) {
        if (gcStatNum[i] > 0) {
            meanGC[i] = meanGC[i] / 255.0 / gcStatNum[i];
        }
    }

    delete[] gcStatNum;

    if (totalNum == 0)
        return 0.0;
    else
        return (double) dupNum / (double) totalNum;
}