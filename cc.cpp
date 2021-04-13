
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>

using namespace std;

int main() {

    string fn[4] = {"preStatsTotalBase", "preStatsTotalQual", "postStatsTotalBase", "postStatsTotalQual"};
    string fr[4] = {"../STD/preStatsTotalBase", "../STD/preStatsTotalQual", "../STD/postStatsTotalBase",
                    "../STD/postStatsTotalQual"};

    fstream out, in, inn;

    auto tlen = 100;
    for (int tt = 0; tt < 4; tt++) {
        uint *f1 = new uint[tlen];
        long *f2 = new long[tlen];

        cout << fn[tt] << " " << fr[tt] << endl;

        in.open(fn[tt].c_str(), ios::in | ios::binary);
        inn.open(fr[tt].c_str(), ios::in | ios::binary);
        if (!in) {
            printf("Can't open file \"%s\"\n", fn[tt].c_str());
        } else if (!inn) {
            printf("Can't open file \"%s\"\n", fr[tt].c_str());
        } else {
            in.seekg(0, ios::beg);
            in.read(reinterpret_cast<char *>(f1), tlen * sizeof(uint));
            inn.seekg(0, ios::beg);
            inn.read(reinterpret_cast<char *>(f2), tlen * sizeof(long));
            printf("=================================================\n");
            for (int i = 0; i < tlen; i++) {
                if (f1[i] != f2[i]) {
                    printf("GG on test %d\n", i);
                    break;
                }
            }
            printf("=================================================\n");
        }
        in.close();
        inn.close();
    }
    return 0;
}