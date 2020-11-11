#include "basecorrector.h"
#include "util.h"

BaseCorrector::BaseCorrector(){
}


BaseCorrector::~BaseCorrector(){
}

int BaseCorrector::correctByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr) {
    OverlapResult ov = OverlapAnalysis::analyze(r1, r2);
    return correctByOverlapAnalysis(r1, r2, fr, ov);
}

int BaseCorrector::correctByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr, OverlapResult ov) {
    // we only correct overlap with up to 5 mismatches
    if(ov.diff == 0 || ov.diff > 5)
        return 0;

    int ol = ov.overlap_len;
    int start1 = max(0, ov.offset);
    int start2 = r2->length() -  max(0, -ov.offset) - 1;

    const char* seq1 = r1->mSeq.c_str();
    const char* seq2 = r2->mSeq.c_str();
    const char* qual1 = r1->mQuality.c_str();
    const char* qual2 = r2->mQuality.c_str();

    const char GOOD_QUAL = num2qual(30);
    const char BAD_QUAL = num2qual(14);

    int corrected = 0;
    int uncorrected = 0;
    bool r1Corrected = false;
    bool r2Corrected = false;
    for(int i=0; i<ol; i++) {
        int p1 = start1 + i;
        int p2 = start2 - i;

        if(seq1[p1] != complement(seq2[p2])) {
            if(qual1[p1] >= GOOD_QUAL && qual2[p2] <= BAD_QUAL) {
                // use R1
                //TODO WRONG
                r2->mSeq.mstr[p2+r2->mSeq.front]=complement(seq1[p1]);
                r2->mQuality.mstr[p2+r2->mQuality.len] = qual1[p1];
                corrected++;
                r2Corrected = true;
                if(fr) {
                    fr->addCorrection(seq2[p2], complement(seq1[p1]));
                }
            } else if(qual2[p2] >= GOOD_QUAL && qual1[p1] <= BAD_QUAL) {
                // use R2
                //TODO WRONG
                r1->mSeq.mstr[p1+r1->mSeq.front]=complement(seq2[p2]);
                r1->mQuality.mstr[p1+r1->mQuality.len] = qual2[p2];
                corrected++;
                r1Corrected = true;
                if(fr) {
                    fr->addCorrection(seq1[p1], complement(seq2[p2]));
                }
            } else {
                uncorrected++;
            }
        }
    }

    // should never happen
    if(uncorrected + corrected != ov.diff) {
        static bool warned = false;
        if(!warned){
            cerr << "WARNING: the algorithm is wrong! uncorrected + corrected != ov.diff" << endl;
            warned = true;
        }
    }

    if(corrected > 0 && fr) {
        if(r1Corrected && r2Corrected)
            fr->incCorrectedReads(2);
        else
            fr->incCorrectedReads(1);
    }

    return corrected;
}

bool BaseCorrector::test() {
    Read r1("@name",
        "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCACGGGG",
        "+",
        "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/EEEEE");
    Read r2("@name",
        "AAAAAAAAAACCCCGGGGAAAATTTTAAAATTGGGGGGGGGGTGGGGGGGGGGGGG",
        "+",
        "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/EEEEEEEEEEEEE");

    correctByOverlapAnalysis(&r1, &r2, NULL);

    if(r1.mSeq != "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG")
        return false;
    if(r2.mSeq != "AAAAAAAAAACCCCGGGGAAAATTTTAAAATTGGGGGGGGGGGGGGGGGGGGGGGG")
        return false;
    if(r1.mQuality != "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE")
        return false;
    if(r2.mQuality != "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE")
        return false;

    return true;
}