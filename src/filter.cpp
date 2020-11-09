#include "processor.h"
#include "peprocessor.h"
#include "seprocessor.h"
#include "overlapanalysis.h"

Filter::Filter(Options* opt){
    mOptions = opt;
}


Filter::~Filter(){
}

int Filter::passFilter(Read* r) {
    if(r == NULL || r->length()==0) {
        return FAIL_LENGTH;
    }

    int rlen = r->length();
    int lowQualNum = 0;
    int nBaseNum = 0;

    // need to recalculate lowQualNum and nBaseNum if the corresponding filters are enabled
    if(mOptions->qualfilter.enabled || mOptions->lengthFilter.enabled) {
        const char* seqstr = r->mSeq.mStr.c_str();
        const char* qualstr = r->mQuality.c_str();

        for(int i=0; i<rlen; i++) {
            char base = seqstr[i];
            char qual = qualstr[i];

            if(qual < mOptions->qualfilter.qualifiedQual)
                lowQualNum ++;

            if(base == 'N')
                nBaseNum++;
        }
    }

    if(mOptions->qualfilter.enabled) {
        if(lowQualNum > (mOptions->qualfilter.unqualifiedPercentLimit * rlen / 100.0) )
            return FAIL_QUALITY;
        else if(nBaseNum > mOptions->qualfilter.nBaseLimit )
            return FAIL_N_BASE;
    }

    if(mOptions->lengthFilter.enabled) {
        if(rlen < mOptions->lengthFilter.requiredLength)
            return FAIL_LENGTH;
        if(mOptions->lengthFilter.maxLength > 0 && rlen > mOptions->lengthFilter.maxLength)
            return FAIL_TOO_LONG;
    }

    if(mOptions->complexityFilter.enabled) {
        if(!passLowComplexityFilter(r))
            return FAIL_COMPLEXITY;
    }

    return PASS_FILTER;
}

bool Filter::passLowComplexityFilter(Read* r) {
    int diff = 0;
    int length = r->length();
    if(length <= 1)
        return false;
    const char* data = r->mSeq.mStr.c_str();
    for(int i=0; i<length-1; i++) {
        if(data[i] != data[i+1])
            diff++;
    }
    if( (double)diff/(double)(length-1) >= mOptions->complexityFilter.threshold )
        return true;
    else
        return false;
}

Read* Filter::trimAndCut(Read* r, int front, int tail) {
    // return the same read for speed if no change needed
    if(front == 0 && tail == 0 && !mOptions->qualityCut.enabled5 && !mOptions->qualityCut.enabled3)
        return r;


    int rlen = r->length() - front - tail ; 
    if (rlen < 0)
        return NULL;

    if(front == 0 && !mOptions->qualityCut.enabled5 && !mOptions->qualityCut.enabled3){
        r->resize(rlen);
        return r;
    } else if(!mOptions->qualityCut.enabled5 && !mOptions->qualityCut.enabled3){
        r->mSeq.mStr._substr(front, rlen);
        r->mQuality = r->mQuality.substr(front, rlen);
        return r;
    }

    // need quality cutting

    int w = mOptions->qualityCut.windowSize;
    int l = r->length();
    const char* qualstr = r->mQuality.c_str();
    const char* seq = r->mSeq.mStr.c_str();
    // quality cutting forward
    if(mOptions->qualityCut.enabled5) {
        int s = front;
        if(l - front - tail - w <= 0)
            return NULL;

        int totalQual = 0;

        // preparing rolling
        for(int i=0; i<w-1; i++)
            totalQual += qualstr[s+i];

        for(s=front; s+w<l-tail; s++) {
            totalQual += qualstr[s+w-1];
            // rolling
            if(s > front) {
                totalQual -= qualstr[s-1];
            }
            // add 33 for phred33 transforming
            if((double)totalQual / (double)w >= 33 + mOptions->qualityCut.quality)
                break;
        }

        // the trimming in front is forwarded and rlen is recalculated
        if(s >0 )
            s = s+w-1;
        while(s<l && seq[s] == 'N')
            s++;
        front = s;
        rlen = l - front - tail;
    }

    // quality cutting backward
    if(mOptions->qualityCut.enabled3) {
        if(l - front - tail - w <= 0)
            return NULL;

        int totalQual = 0;
        int t = l - tail - 1;

        // preparing rolling
        for(int i=0; i<w-1; i++)
            totalQual += qualstr[t-i];

        for(t=l-tail-1; t-w>=front; t--) {
            totalQual += qualstr[t-w+1];
            // rolling
            if(t < l-tail-1) {
                totalQual -= qualstr[t+1];
            }
            // add 33 for phred33 transforming
            if((double)totalQual / (double)w >= 33 + mOptions->qualityCut.quality)
                break;
        }

        if(t < l-1)
            t = t-w+1;
        while(t>=0 && seq[t] == 'N')
            t--;
        rlen = t - front + 1;
    }

    if(rlen <= 0 || front >= l-1)
        return NULL;

    r->mSeq.mStr._substr(front, rlen);
    r->mQuality = r->mQuality.substr(front, rlen);

    return r;
}

bool Filter::filterByIndex(Read* r) {
    if(mOptions->indexFilter.enabled) {
        if( match(mOptions->indexFilter.blacklist1, r->firstIndex(), mOptions->indexFilter.threshold) )
            return true;
    }
    return false;
}

bool Filter::filterByIndex(Read* r1, Read* r2) {
    if(mOptions->indexFilter.enabled) {
        if( match(mOptions->indexFilter.blacklist1, r1->firstIndex(), mOptions->indexFilter.threshold) )
            return true;
        if( match(mOptions->indexFilter.blacklist2, r2->lastIndex(), mOptions->indexFilter.threshold) )
            return true;
    }
    return false;
}

bool Filter::match(vector<string>& list, string target, int threshold) {
    for(int i=0; i<list.size(); i++) {
        int diff = 0;
        int len1 = list[i].length();
        int len2 = target.length();
        for(int s=0; s<len1 && s<len2; s++) {
            if(list[i][s] != target[s]) {
                diff++;
                if(diff>threshold)
                    break;
            }
        }
        if(diff <= threshold)
            return true;
    }
    return false;
}

bool Filter::test() {
    Read r("@name",
        "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTT",
        "+",
        "/////CCCCCCCCCCCC////CCCCCCCCCCCCCC////E");
    Options opt;
    opt.qualityCut.enabled5 = true;
    opt.qualityCut.enabled3 = true;
    opt.qualityCut.windowSize = 4;
    opt.qualityCut.quality = 20;
    Filter filter(&opt);
    Read* ret = filter.trimAndCut(&r, 0, 1);
    ret->print();
    
    return ret->mSeq.mStr == "CCCCCCCCCCCCCCCCCCCCCCCCCCCC"
        && ret->mQuality == "CCCCCCCCCCC////CCCCCCCCCCCCC";
}