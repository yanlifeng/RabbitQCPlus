#include "overlapanalysis.h"

#ifdef Vec512

#include<immintrin.h>

#endif

#ifdef Vec256

#include<immintrin.h>

#endif

OverlapAnalysis::OverlapAnalysis() {
}


OverlapAnalysis::~OverlapAnalysis() {
}

OverlapResult OverlapAnalysis::analyze(Read *r1, Read *r2, int overlapDiffLimit, int overlapRequire) {
    return analyze(r1->mSeq, r2->mSeq, overlapDiffLimit, overlapRequire);
}


#ifdef Vec512

// ported from the python code of AfterQC
OverlapResult OverlapAnalysis::analyze(Sequence &r1, Sequence &r2, int overlapDiffLimit, int overlapRequire) {
    Sequence rcr2 = ~r2;
    int len1 = r1.length();
    int len2 = rcr2.length();
    // use the pointer directly for speed
    const char *str1 = r1.mStr.c_str();
    const char *str2 = rcr2.mStr.c_str();

    int complete_compare_require = 50;

    int overlap_len = 0;
    int offset = 0;
    int diff = 0;

    // forward
    // a match of less than overlapRequire is considered as unconfident
    while (offset < len1 - overlapRequire) {
        // the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1 - offset, len2);
        diff = 0;
        int leng = min(complete_compare_require, overlap_len);
        for (int i = 0; i < leng; i += 64) {
            uint64 tag = (1ll << min(64, leng - i)) - 1;
            __m512i t1 = _mm512_maskz_loadu_epi8(tag, str1 + offset + i);
            __m512i t2 = _mm512_maskz_loadu_epi8(tag, str2 + i);
            __mmask64 res = _mm512_cmp_epi8_mask(t1, t2, 4);
            diff += _mm_popcnt_u64(res);
            if (diff > overlapDiffLimit)break;
        }

        if (diff <= overlapDiffLimit) {
            OverlapResult ov;
            ov.overlapped = true;
            ov.offset = offset;
            ov.overlap_len = overlap_len;
            ov.diff = diff;
            return ov;
        }
        offset += 1;
    }


    // reverse
    // in this case, the adapter is sequenced since TEMPLATE_LEN < SEQ_LEN
    // check if distance can get smaller if offset goes negative
    // this only happens when insert DNA is shorter than sequencing read length, and some adapter/primer is sequenced but not trimmed cleanly
    // we go reversely
    offset = 0;
    while (offset > -(len2 - overlapRequire)) {
        // the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1, len2 - abs(offset));

        diff = 0;

        int leng = min(complete_compare_require, overlap_len);
        for (int i = 0; i < leng; i += 64) {
            uint64 tag = (1ll << min(64, leng - i)) - 1;
            __m512i t1 = _mm512_maskz_loadu_epi8(tag, str1 + i);
            __m512i t2 = _mm512_maskz_loadu_epi8(tag, str2 - offset + i);
            __mmask64 res = _mm512_cmp_epi8_mask(t1, t2, 4);
            diff += _mm_popcnt_u64(res);
            if (diff > overlapDiffLimit)break;
        }
        if (diff <= overlapDiffLimit) {
            OverlapResult ov;
            ov.overlapped = true;
            ov.offset = offset;
            ov.overlap_len = overlap_len;
            ov.diff = diff;
            return ov;
        }
        offset -= 1;
    }

    OverlapResult ov;
    ov.overlapped = false;
    ov.offset = ov.overlap_len = ov.diff = 0;
    return ov;
}

#elif Vec256

OverlapResult OverlapAnalysis::analyze(Sequence &r1, Sequence &r2, int overlapDiffLimit, int overlapRequire) {
    Sequence rcr2 = ~r2;
    int len1 = r1.length();
    int len2 = rcr2.length();
    // use the pointer directly for speed
    const char *str1 = r1.mStr.c_str();
    const char *str2 = rcr2.mStr.c_str();

    int complete_compare_require = 50;

    int overlap_len = 0;
    int offset = 0;
    int diff = 0;

    // forward
    // a match of less than overlapRequire is considered as unconfident
    while (offset < len1 - overlapRequire) {
        // the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1 - offset, len2);
        diff = 0;
        int leng = min(complete_compare_require, overlap_len);
        int i = 0;
        for (; i + 32 <= leng; i += 32) {
            __m256i t1 = _mm256_loadu_si256(reinterpret_cast<const __m256i_u *>(str1 + offset + i));
            __m256i t2 = _mm256_loadu_si256(reinterpret_cast<const __m256i_u *>(str2 + i));
            __m256i res = _mm256_cmpeq_epi8(t1, t2);
            for (int j = 0; j < 32; j++) {
                if (res[j] == 0)diff++;
            }
            if (diff > overlapDiffLimit)break;
        }
        for (; i < leng; i++) {
            if (str1[offset + i] != str2[i]) {
                diff += 1;
                if (diff > overlapDiffLimit)break;
            }
        }

        if (diff <= overlapDiffLimit) {
            OverlapResult ov;
            ov.overlapped = true;
            ov.offset = offset;
            ov.overlap_len = overlap_len;
            ov.diff = diff;
            return ov;
        }
        offset += 1;
    }


    // reverse
    // in this case, the adapter is sequenced since TEMPLATE_LEN < SEQ_LEN
    // check if distance can get smaller if offset goes negative
    // this only happens when insert DNA is shorter than sequencing read length, and some adapter/primer is sequenced but not trimmed cleanly
    // we go reversely
    offset = 0;
    while (offset > -(len2 - overlapRequire)) {
        // the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1, len2 - abs(offset));

        diff = 0;
        int leng = min(complete_compare_require, overlap_len);


        int i = 0;
        for (; i + 32 <= leng; i += 32) {
            __m256i t1 = _mm256_loadu_si256(reinterpret_cast<const __m256i_u *>(str1 + i));
            __m256i t2 = _mm256_loadu_si256(reinterpret_cast<const __m256i_u *>(str2 - offset + i));
            __m256i res = _mm256_cmpeq_epi8(t1, t2);
            for (int j = 0; j < 32; j++) {
                if (res[j] == 0)diff++;
            }
            if (diff > overlapDiffLimit)break;
        }
        for (; i < leng; i++) {
            if (str1[i] != str2[-offset + i]) {
                diff += 1;
                if (diff > overlapDiffLimit)break;
            }
        }
        if (diff <= overlapDiffLimit) {
            OverlapResult ov;
            ov.overlapped = true;
            ov.offset = offset;
            ov.overlap_len = overlap_len;
            ov.diff = diff;
            return ov;
        }
        offset -= 1;
    }

    OverlapResult ov;
    ov.overlapped = false;
    ov.offset = ov.overlap_len = ov.diff = 0;
    return ov;
}

#else
static char reMap[123] = {'0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
                   '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
                   '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
                   '0', '0', '0', '0', '0', 'T', 'B', 'G', 'D', 'E', 'F', 'C', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
                   'P', 'Q', 'R', 'S', 'A', 'U', 'V', 'W', 'X', 'Y', 'Z', '0', '0', '0', '0', '0', '0', 'T', 'b', 'G',
                   'd', 'e', 'f', 'C', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 'A', 'u', 'v', 'w',
                   'x', 'y', 'z'};

// ported from the python code of AfterQC
OverlapResult OverlapAnalysis::analyze(Sequence &r1, Sequence &r2, int overlapDiffLimit, int overlapRequire) {
    int len1 = r1.length();
    int len2 = r2.length();
    // use the pointer directly for speed
    const char *str1 = r1.mStr.c_str();
    const char *str2 = r2.mStr.c_str();


    int complete_compare_require = 50;

    int overlap_len = 0;
    int offset = 0;
    int diff = 0;

    // forward
    // a match of less than overlapRequire is considered as unconfident
    while (offset < len1 - overlapRequire) {
        // the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1 - offset, len2);

        diff = 0;
        int i = 0;
        for (; i < overlap_len; i++) {
            if (str1[offset + i] != reMap[str2[len2 - 1 - i]]) {
                diff += 1;
                if (diff == overlapDiffLimit && i < complete_compare_require)
                    break;
            }
        }

        if (diff < overlapDiffLimit || (diff >= overlapDiffLimit && i >= complete_compare_require)) {
            OverlapResult ov;
            ov.overlapped = true;
            ov.offset = offset;
            ov.overlap_len = overlap_len;
            ov.diff = diff;
            return ov;
        }

        offset += 1;
    }


    // reverse
    // in this case, the adapter is sequenced since TEMPLATE_LEN < SEQ_LEN
    // check if distance can get smaller if offset goes negative
    // this only happens when insert DNA is shorter than sequencing read length, and some adapter/primer is sequenced but not trimmed cleanly
    // we go reversely
    offset = 0;
    while (offset > -(len2 - overlapRequire)) {
        // the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1, len2 - abs(offset));

        diff = 0;
        int i = 0;
        for (i = 0; i < overlap_len; i++) {
            if (str1[i] != reMap[str2[len2 - 1 + offset - i]]) {
                diff += 1;
                if (diff >= overlapDiffLimit && i < complete_compare_require)
                    break;
            }
        }

        if (diff < overlapDiffLimit || (diff >= overlapDiffLimit && i >= complete_compare_require)) {
            OverlapResult ov;
            ov.overlapped = true;
            ov.offset = offset;
            ov.overlap_len = overlap_len;
            ov.diff = diff;
            return ov;
        }

        offset -= 1;
    }

    OverlapResult ov;
    ov.overlapped = false;
    ov.offset = ov.overlap_len = ov.diff = 0;
    return ov;
}

#endif

bool OverlapAnalysis::test() {
    //Sequence r1("CAGCGCCTACGGGCCCCTTTTTCTGCGCGACCGCGTGGCTGTGGGCGCGGATGCCTTTGAGCGCGGTGACTTCTCACTGCGTATCGAGCCGCTGGAGGTCTCCC");
    //Sequence r2("ACCTCCAGCGGCTCGATACGCAGTGAGAAGTCACCGCGCTCAAAGGCATCCGCGCCCACAGCCACGCGGTCGCGCAGAAAAAGGGGCCCGTAGGCGCGGCTCCC");

    Sequence r1("CAGCGCCTACGGGCCCCTTTTTCTGCGCGACCGCGTGGCTGTGGGCGCGGATGCCTTTGAGCGCGGTGACTTCTCACTGCGTATCGAGC");
    Sequence r2("ACCTCCAGCGGCTCGATACGCAGTGAGAAGTCACCGCGCTCAAAGGCATCCGCGCCCACAGCCACGCGGTCGCGCAGAAAAAGGGGTCC");

    OverlapResult ov = OverlapAnalysis::analyze(r1, r2);

    return ov.overlapped && ov.offset == 10 && ov.overlap_len == 79 && ov.diff == 1;
}