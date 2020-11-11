#ifndef STRINGPOOL_H
#define STRINGPOOL_H
#include <cstdio>
#include <chrono>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include <stack>
#include <set>
#include "util.h"

/*
 *  出现修改mstr的操作的时候，一定需要添加front来进行偏移
 */
class StringPool{
public:
    StringPool();
    StringPool(std::string s);


    char operator [](int i);


    std::string operator +   (const std::string &b);
    std::string operator =   (const std::string &b);
    bool        operator ==  (const std::string &b);
    bool        operator !=  (const std::string &b);
    bool        operator ==  (StringPool b);
    bool        operator !=  (StringPool b);
    StringPool  operator ~   ();

    //用于更换char
    void set(int pos,char x);


    StringPool reverseComplement();
    StringPool reverse();

    // 生成新的substring
    std::string substr(int pos);
    std::string substr(int pos,int n);


    // 在原始的字符串上进行操作，不生成新的substring
    void _substr(int pos);

    //resize 看着是read类上已经进行了判断，这里还是再判断一次吧
    void resize(int n);


    void _substr(int pos,int n);

    std::string toString();
    char* c_str();
    int length();
    int size();

public:
    std::string mstr; // 字符串(可能出现赘余)
    int front; // 字符串开始的位置
    int len; // 字符串的长度
};
#endif
