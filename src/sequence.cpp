#include "sequence.h"

Sequence::Sequence() {
}

Sequence::Sequence(string seq) {
    mStr = seq;
}

void Sequence::print() {
    std::cerr << mStr;
}

int Sequence::length() {
    return mStr.length();
}

static char reMap[123] = {'0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
                          '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
                          '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
                          '0', '0', '0', '0', '0', '0', '0', '0', 'T', 'B', 'G', 'D', 'E', 'F', 'C', 'H', 'I', 'J', 'K',
                          'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'A', 'U', 'V', 'W', 'X', 'Y', 'Z', '0', '0', '0', '0',
                          '0', '0', 'T', 'b', 'G', 'd', 'e', 'f', 'C', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q',
                          'r', 's', 'A', 'u', 'v', 'w', 'x', 'y', 'z'};

Sequence Sequence::reverseComplement() {
    string str(mStr.length(), 0);
    for (int c = 0; c < mStr.length(); c++) {
        str[mStr.length() - c - 1] = reMap[mStr[c]];
    }
    return Sequence(str);
}

Sequence Sequence::operator~() {
    return reverseComplement();
}

bool Sequence::test() {
    Sequence s("AAAATTTTCCCCGGGG");
    Sequence rc = ~s;
    if (s.mStr != "AAAATTTTCCCCGGGG") {
        cerr << "Failed in reverseComplement() expect AAAATTTTCCCCGGGG, but get " << s.mStr;
        return false;
    }
    if (rc.mStr != "CCCCGGGGAAAATTTT") {
        cerr << "Failed in reverseComplement() expect CCCCGGGGAAAATTTT, but get " << rc.mStr;
        return false;
    }
    return true;
}