#include "StringPool.h"

StringPool::StringPool(){};

StringPool::StringPool(std::string s){
    this->mstr = s;
    this->front = 0;
    this->len = s.size();
}


char StringPool::operator [](int i){
    if (i>=this->len){return this->mstr[this->front];}
    return this->mstr[this->front+i];
}


std::string StringPool::operator + (const std::string &b){
    return std::string(this->mstr.substr(this->front,this->len)+b);
}


bool StringPool::operator !=(const std::string &b){
    bool fg=false;
    if (this->len != b.size()) return true;
    for (int i=0;i<this->len;i++) if (b[i]!=this->mstr[i+this->front]){fg=true;break;}
    return fg;
}

bool StringPool::operator ==(const std::string &b){
    bool fg=true;
    if (this->len != b.size()) return false;
    for (int i=0;i<this->len;i++) if (b[i]!=this->mstr[i+this->front]){fg=false;break;}
    return fg;
}

// 生成新的substd::string
std::string StringPool::substr(int pos){
    if (pos >= this->len || pos<0) return this->mstr.substr(this->front,this->len);
    return this->mstr.substr(pos+this->front,this->len-pos);
}
std::string StringPool::substr(int pos,int n){
    if ( pos+n >= this->len || pos<0) return this->mstr.substr(this->front,this->len);
    return this->mstr.substr(pos + this->front,n);
}


// 在原始的字符串上进行操作，不生成新的substring
void StringPool::_substr(int pos){
    if (pos >= this->len || pos<0) return;
    this->len -= pos;
    this->front += pos;
}

//resize 看着是read类上已经进行了判断，这里还是再判断一次吧
void StringPool::resize(int n){
    if (this->len <=n || n<0) return;
    this->len=n;
}


void StringPool::_substr(int pos,int n){
    if ( pos+n >= this->len || pos<0 ) return;;
    this->len = n;
    this->front += pos;
}

std::string StringPool::toString(){
    return this->mstr.substr(this->front,this->len);
}
char* StringPool::c_str(){return &this->mstr[front];}
int StringPool::length(){return this->len;}
//int main(int argc,char **argv){
//    StringPool s("AAAATTTTCCCCGGGG");
//    if (s != "AAAATTTTCCCGGGG"){
//        cout << "okkkk" << endl;
//    }
//}

