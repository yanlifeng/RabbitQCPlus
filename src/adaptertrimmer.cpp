#include "adaptertrimmer.h"

#ifdef Vec512

#include<immintrin.h>

#endif

AdapterTrimmer::AdapterTrimmer() {
}


AdapterTrimmer::~AdapterTrimmer() {
}

bool AdapterTrimmer::trimByOverlapAnalysis(Read *r1, Read *r2, FilterResult *fr) {
    OverlapResult ov = OverlapAnalysis::analyze(r1, r2);
    return trimByOverlapAnalysis(r1, r2, fr, ov);
}

bool AdapterTrimmer::trimByOverlapAnalysis(Read *r1, Read *r2, FilterResult *fr, OverlapResult ov) {
    int ol = ov.overlap_len;
    if (ov.diff <= 5 && ov.overlapped && ov.offset < 0 && ol > r1->length() / 3) {
        string adapter1 = r1->mSeq.mStr.substr(ol, r1->length() - ol);
        string adapter2 = r2->mSeq.mStr.substr(ol, r2->length() - ol);

        if (_DEBUG) {
            cerr << adapter1 << endl;
            cerr << adapter2 << endl;
            cerr << "overlap:" << ov.offset << "," << ov.overlap_len << ", " << ov.diff << endl;
            r1->print();
            r2->reverseComplement()->print();
            cerr << endl;
        }

        r1->mSeq.mStr = r1->mSeq.mStr.substr(0, ol);
        r1->mQuality = r1->mQuality.substr(0, ol);
        r2->mSeq.mStr = r2->mSeq.mStr.substr(0, ol);
        r2->mQuality = r2->mQuality.substr(0, ol);

        fr->addAdapterTrimmed(adapter1, adapter2);
        return true;
    }
    return false;
}

bool AdapterTrimmer::trimBySequence(Read *r, FilterResult *fr, string &adapterseq, bool isR2) {
    const int matchReq = 4;
    const int allowOneMismatchForEach = 8;

    int rlen = r->length();
    int alen = adapterseq.length();

    const char *adata = adapterseq.c_str();
    const char *rdata = r->mSeq.mStr.c_str();

    if (alen < matchReq)
        return false;

    int pos = 0;
    bool found = false;
    int start = 0;
    if (alen >= 16)
        start = -4;
    else if (alen >= 12)
        start = -3;
    else if (alen >= 8)
        start = -2;
    // we start from negative numbers since the Illumina adapter dimer usually have the first A skipped as A-tailing


#ifdef Vec512
    for (pos = start; pos < rlen - matchReq; pos++) {
        int cmplen = min(rlen - pos, alen);
        int allowedMismatch = cmplen / allowOneMismatchForEach;
        int mismatch = 0;
        int l = max(0, -pos);
        for (int i = l; i < cmplen; i += 64) {
            uint64 tag = (1ll << min(64, cmplen - i)) - 1;
            __m512i t1 = _mm512_maskz_loadu_epi8(tag, adata + i);
            __m512i t2 = _mm512_maskz_loadu_epi8(tag, rdata + i + pos);
            __mmask64 res = _mm512_cmp_epi8_mask(t1, t2, 4);
            mismatch += _mm_popcnt_u64(res);
            if (mismatch > allowedMismatch)break;
        }
//
//        for (int i = l; i < cmplen; i += 32) {
//            uint32 tag = (1ll << min(32, cmplen - i)) - 1;
//            __m256i t1 = _mm256_maskz_loadu_epi8(tag, adata + i);
//            __m256i t2 = _mm256_maskz_loadu_epi8(tag, rdata + i + pos);
//            __mmask32 res = _mm256_cmp_epi8_mask(t1, t2, 4);
//            mismatch += _mm_popcnt_u32(res);
//            if (mismatch > allowedMismatch)break;
//        }

//        for (int i = l; i < cmplen; i += 16) {
//            uint16 tag = (1 << min(16, cmplen - i)) - 1;
//            __m128i t1 = _mm_maskz_loadu_epi8(tag, adata + i);
//            __m128i t2 = _mm_maskz_loadu_epi8(tag, rdata + i + pos);
//            __mmask16 res = _mm_cmp_epi8_mask(t1, t2, 4);
//            mismatch += _mm_popcnt_u32(res);
//            if (mismatch > allowedMismatch)break;
//        }

        if (mismatch <= allowedMismatch) {
            found = true;
            break;
        }
    }
#else
    for (pos = start; pos < rlen - matchReq; pos++) {
        int cmplen = min(rlen - pos, alen);
        int allowedMismatch = cmplen / allowOneMismatchForEach;
        int mismatch = 0;
        bool matched = true;
        for (int i = max(0, -pos); i < cmplen; i++) {
            mismatch += adata[i] != rdata[i + pos];
            if (mismatch > allowedMismatch)break;
        }
        if (mismatch <= allowedMismatch) {
            found = true;
            break;
        }
    }
#endif
//TODO
    if (found) {
        if (pos < 0) {
            string adapter = adapterseq.substr(0, alen + pos);
            r->mSeq.mStr.resize(0);
            r->mQuality.resize(0);
            if (fr) {
                fr->addAdapterTrimmed(adapter, isR2);
            }
        } else {
            string adapter = r->mSeq.mStr.substr(pos, rlen - pos);
            r->mSeq.mStr = r->mSeq.mStr.substr(0, pos);
            r->mQuality = r->mQuality.substr(0, pos);
            if (fr) {
                fr->addAdapterTrimmed(adapter, isR2);
            }
        }
        return true;
    }

    return false;
}

bool AdapterTrimmer::test() {
    Read r("@name",
           "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG",
           "+",
           "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////E");
    string adapter = "TTTTCCACGGGGATACTACTG";
    bool trimmed = AdapterTrimmer::trimBySequence(&r, NULL, adapter);
    return r.mSeq.mStr == "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAA";
}