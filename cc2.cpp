
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>

using namespace std;

int main() {

    string fn[8] = {"preStats1TotalBase", "preStats2TotalBase", "preStats1TotalQual", "preStats2TotalQual",
                    "postStats1TotalBase", "postStats2TotalBase", "postStats1TotalQual", "postStats2TotalQual"};
    string fr[8] = {"../STD/preStats1TotalBase", "../STD/preStats2TotalBase", "../STD/preStats1TotalQual",
                    "../STD/preStats2TotalQual", "../STD/postStats1TotalBase", "../STD/postStats2TotalBase",
                    "../STD/postStats1TotalQual", "../STD/postStats2TotalQual"};

    fstream out, in, inn;

    auto tlen = 100;
    for (int tt = 0; tt < 8; tt++) {
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