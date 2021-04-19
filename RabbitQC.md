

## RabbitQC

TODO
  - [ ] 研究main函数中每一个参数的意义
  - [x] 打开trimAndCut
  - [x] 打开polg
  - [x]  解决mDuplicate线程不安全问题
  - [x]  打开mUmiProcessor
  - [ ]  --[haoz:] I think it 3GS can add here
  - [ ]  解决带adapter的数据每次统计信息输出都不一样的问题
  - [x]  check all！
  - [ ]  把string换成char*  不在nwe新的空间
  - [ ]  关闭多线程动态调度，静态分配然后多线程处理，check多线程正确性
  - [x]  使用带adapter的数据
  - [x]  去掉out_str
  - [ ]  改造Read
  - [x] 关于空间
  - [ ] 解决报错./rabbit_qc -w 1 -i ../data/SRR2496709_1.fastq -o p.fq -t 4 -5
  - [x] 为啥运行空间稳定1g不变，
  - [ ] 优化duplicate中的round
  - [ ] state中尝试*8/16的向量化
  - [ ]  


#### 1105

开始动kmer了。

把base2val函数去掉了，换成一个静态的数组

|                                 | total | Stats::statRead |
| ------------------------------- | ----- | --------------- |
| Init                            | 10.78 | 2.43            |
| 把统计和kmer分开&&取消了val函数 | 10.63 | 2.32            |
|                                 |       |                 |

#### 1106

摸大鱼

#### 1107

原来的check方式不太可，因为数据没有被过滤掉的，所以不论怎么改kmer，输入和输出都是一样的。

修改之后的check方法是把kmer计算时几个关键的数组输出一下，然后check这个数组是否不变。

#### 1108

今天为了statRead函数里面的向量化做了一些准备工作，把for去除掉了。

但是发现statRead的两个循环都不适合做向量化，就仅仅简化了一下代码；速度上的话，对于N很多的序列有很好的优化效果，但是似乎很少有N很多的序列，所以影响不大。

|                                                             | total | Stats::statRead |
| ----------------------------------------------------------- | ----- | --------------- |
| Init                                                        | 10.78 | 2.43            |
| 把统计和kmer分开&&取消了val函数                             | 10.63 | 2.32            |
| 又合起来了&&修改了kmer计算函数+-funroll-loops -flto  -mavx2 | 9.8   | 2.13            |

#### 1109

把代码整理了一下，学了学git传上去了。

开始看Duplicate::statRead了，把Duplicate::seq2int里面的switch优化掉了。

|                                                             | total | Stats::statRead | Duplicate::statRead |
| ----------------------------------------------------------- | ----- | --------------- | ------------------- |
| Init                                                        | 10.78 | 2.43            | 1.79                |
| 把统计和kmer分开&&取消了val函数                             | 10.63 | 2.32            | 1.79                |
| 又合起来了&&修改了kmer计算函数+-funroll-loops -flto  -mavx2 | 9.80  | 2.13            | 1.79                |
| 把Duplicate::seq2int里面的switch优化掉                      | 9.09  | 2.04            | 0.60                |

#### 1110 1111 1112

入党材料治我

#### 1113

换了些数据重新跑

##### Init

| SRR2496709_1 | SRR2496709_2 | SRR2530740 | SRR2530740.sra_1 | SRR2530740.sra_2 |
| ------------ | ------------ | ---------- | ---------------- | ---------------- |
| 3.3G         | 3.3G         | 7.5G       | 6.6G             | 6.6G             |
| Adapter Yes  | Adapter Yes  | Adapter No | Adapter No       | Adapter No       |
| 47           | 49           | 74         | 66               | 67               |
| 52           | 54           | 89         | 80               | 79               |



##### 1113V1

|      |      |      |      |      |
| ---- | ---- | ---- | ---- | ---- |
| 48.5 | 49.9 | 75.5 | 66.1 | 67.0 |
|      |      |      |      |      |
|      |      |      |      |      |

##### 1203

昨天突然发现mDuplicate好像不是线程安全的，为了进一步确认，准备多线程跑几次，看看是不是答案不同。但是因为多线程随机分任务，每次运行的结果本来就不相同。

##### 1204

昨天跑着跑着就下课了，然后睡了一下午。。

现在带adapter的数据跑的结果有问题，但是SRR2530740太大了，切了600M的数据下来做测试。

| V    | statReadCost | processSingleEndCost | totleCost |
| ---- | ------------ | -------------------- | --------- |
| STD  | 1.92         | 5.51                 | 7.87      |
| 1204 | 1.83         | 4.12                 | 6.22      |

现在的加速主要来自mDuplicate的statRead，Status里面的statRead优化不大，后面主要从字符串的角度进行优化。

##### 1209

```c++

bool SingleEndProcessor::processSingleEnd(ReadPack *pack, ThreadConfig *config) {


    //debug for losing last line
    //cerr << pack->data[pack->count - 1]->mName << endl;
    //cerr << pack->data[pack->count - 1]->mQuality << endl;
    //debug
    string outstr;
    int readPassed = 0;
    //------------------my thinking---------------------------------
    /*
    if (mOptions -> thirdgene){
      for(int p = 0; p < pack->data[p]; p++){
        Read* or1 = pack->data[p];
        config -> getTGStats() -> tgsStatRead(or1);
      }
    }
    */
    //1. add TGStats class in project
    //2. add getTGStats function in ThreadConfig file;
    //3. add TGStats menber variable in ThreadConfig class
    //---------------------------------------------------
//    printf("Totel : \n");
//    printf("pack->count : %d\n", pack->count);
//    printf("total cost1 : %.5f\n", cost1);
//    printf("total cost2 : %.5f\n", cost2);
//    printf("total cost3 : %.5f\n", cost3);
//    printf("total cost4 : %.5f\n", cost4);
    totCnt += 1;
    double t0, t1;


    t0 = get_wall_time();

    for (int p = 0; p < pack->count; p++) {

        // original read1
        Read *or1 = pack->data[p];
        // stats the original read before trimming
        t1 = get_wall_time();
        config->getPreStats1()->statRead(or1);//cost 12/57
        cost1 += get_wall_time() - t1;


        // handling the duplication profiling
        t1 = get_wall_time();

        //TODO maybe mDuplicate is not thread safe
        if (mDuplicate)
            mDuplicate->statRead(or1);//cost 14/57
        cost2 += get_wall_time() - t1;


        t1 = get_wall_time();
        // filter by index
        if (mOptions->indexFilter.enabled && mFilter->filterByIndex(or1)) {
//            printf("mOptions->indexFilter ...");
            delete or1;
            continue;
        }
        cost3 += get_wall_time() - t1;


        t1 = get_wall_time();
        // umi processing

        //TODO maybe this can be the big hotspot
        if (mOptions->umi.enabled) {
//            printf("mOptions->umi ...");
            mUmiProcessor->process(or1);
        }
        cost4 += get_wall_time() - t1;


        t1 = get_wall_time();
        // trim in head and tail, and apply quality cut in sliding window
        //not in
        Read *r1 = mFilter->trimAndCut(or1, mOptions->trim.front1, mOptions->trim.tail1);
        cost5 += get_wall_time() - t1;


        t1 = get_wall_time();
        if (r1 != NULL) {
            //not in because xx is false
            if (mOptions->polyGTrim.enabled)
                PolyX::trimPolyG(r1, config->getFilterResult(), mOptions->polyGTrim.minLen);
            //not in because xx is false
            if (mOptions->polyXTrim.enabled)
                PolyX::trimPolyX(r1, config->getFilterResult(), mOptions->polyXTrim.minLen);
        }
        cost6 += get_wall_time() - t1;


        t1 = get_wall_time();
        //not in because mOptions->adapter.hasSeqR1 is false
        if (r1 != NULL && mOptions->adapter.enabled && mOptions->adapter.hasSeqR1) {
            AdapterTrimmer::trimBySequence(r1, config->getFilterResult(), mOptions->adapter.sequence);
        }
        cost7 += get_wall_time() - t1;


        t1 = get_wall_time();
        if (r1 != NULL) {
            // not in because mOptions->trim.maxLen1 is 0
            if (mOptions->trim.maxLen1 > 0 && mOptions->trim.maxLen1 < r1->length())
                r1->resize(mOptions->trim.maxLen1);
        }
        cost8 += get_wall_time() - t1;


        t1 = get_wall_time();
        //in ! O(len)
        int result = mFilter->passFilter(r1);
        cost9 += get_wall_time() - t1;


        t1 = get_wall_time();
        config->addFilterResult(result);
        cost10 += get_wall_time() - t1;

        //in !
        if (r1 != NULL && result == PASS_FILTER) {
            t1 = get_wall_time();
            outstr += r1->toString();
            cost11b += get_wall_time() - t1;
            // stats the read after filtering
            t1 = get_wall_time();
            config->getPostStats1()->statRead(r1);//cost 23/57
            cost11 += get_wall_time() - t1;
            readPassed++;
        }


        t1 = get_wall_time();
        delete or1;
        // if no trimming applied, r1 should be identical to or1
        if (r1 != or1 && r1 != NULL)
            delete r1;
        cost12 += get_wall_time() - t1;
//        printf("now %d : \n", p);
//        printf("total cost1 : %.5f\n", cost1);
//        printf("total cost2 : %.5f\n", cost2);
//        printf("total cost3 : %.5f\n", cost3);
//        printf("total cost4 : %.5f\n", cost4);
//        printf("\n");
    }
    cost += get_wall_time() - t0;

    // if splitting output, then no lock is need since different threads write different files
    if (!mOptions->split.enabled)
        mOutputMtx.lock();
    if (mOptions->outputToSTDOUT) {
        fwrite(outstr.c_str(), 1, outstr.length(), stdout);
    } else if (mOptions->split.enabled) {
        // split output by each worker thread
        if (!mOptions->out1.empty())
            config->getWriter1()->writeString(outstr);
    } else {
        if (mLeftWriter) {
            char *ldata = new char[outstr.size()];
            memcpy(ldata, outstr.c_str(), outstr.size());
            mLeftWriter->input(ldata, outstr.size());
        }
    }

    if (!mOptions->split.enabled)
        mOutputMtx.unlock();

    if (mOptions->split.byFileLines)
        config->markProcessed(readPassed);
    else
        config->markProcessed(pack->count);

    //delete pack->data;
    std::vector<Read *>().swap(pack->data);
    delete pack;
    return true;
}
```

这一大段代码就是最主要的处理过程，但是gdb单步发现，很多很多string处理的函数都没有进去，就是功能没有开开

```c++
Filtering result:
reads passed filter: 1979553
reads failed due to low quality: 20180
reads failed due to too many N: 267
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0
```

可以看到只是把质量太低和N太多的数据删掉了，并没有进行N的删除之类的操作。

默认的功能大概就是，比较快速的读进来，做一遍初始的统计，然后做修建和删除：

```c++
        Read *r1 = mFilter->trimAndCut(or1, mOptions->trim.front1, mOptions->trim.tail1);
```

但是，这个功能并没有打开，所以r1和or1是一样的，然后做mFilter->passFilter，测试r1的质量怎么样，如果可以的话，就输出，并记录，所以现在的热点就仅仅是输出和统计。

尝试开启trimAndCut：

```bash
./rabbit_qc -w 1 -i ../data/test.fq -o p.fq -t 4 -5 -3 
```

时间从0.05->0.2，还不是热点，具体的函数定义还要再看看。



#### 1210

解决mDuplicate线程不安全问题：

按理说只要在mDuplicate的变量操作上加上锁就好了，但是怎么check自己加的对不对是个大问题。因为thread unsafe的问题只有在多线程的时候才会出现，但是由于每次相当于消费者随机取生产者准备好的数据，也就是说多线程的答案本来就是随机的。

但是merge之后的统计信息应该是一样的，所以我们不再每个线程都输出自己的统计信息，等merge之后在check。

更新了计时函数（所有的线程一起计时）。





#### 1214

GC的计算有个bug，但是改完之后没什么用处，这个暂时放一下，先试试string换char*的思路。

#### 1215

```bash
STD time use test.fq:
total getPreStats1()->statRead(or1) ====: 1.14592
total mDuplicate->statRead(or1) ========: 1.74557
total mOptions->indexFilter()  =========: 0.04799
total mUmiProcessor->process(or1) ======: 0.15901
total mFilter->trimAndCut() ============: 0.05417
total PolyX::trimPolyG() ===============: 0.15901
total trimBySequence ===================: 0.04800
total r1->resize() =====================: 0.04908
total mFilter->passFilter(r1) ==========: 0.12123
total addFilterResult(result) ==========: 0.05077
total outstr += r1->toString() =========: 0.67996
total getPostStats1()->statRead(r1) ====: 1.07540
total delete r1 ========================: 0.33475
total costTotel ========================: 5.55984
total cost =============================: 6.25348
total  =================================: 186
10.9451 

Plus time use test.fq:
total getPreStats1()->statRead(or1) ====: 1.08299
total mDuplicate->statRead(or1) ========: 0.61184
total mOptions->indexFilter()  =========: 0.04844
total mUmiProcessor->process(or1) ======: 0.15745
total mFilter->trimAndCut() ============: 0.05405
total PolyX::trimPolyG() ===============: 0.15745
total trimBySequence ===================: 0.04884
total r1->resize() =====================: 0.04875
total mFilter->passFilter(r1) ==========: 0.10593
total addFilterResult(result) ==========: 0.05057
total outstr += r1->toString() =========: 0.64141
total getPostStats1()->statRead(r1) ====: 1.03461
total delete r1 ========================: 0.33882
total costTotel ========================: 4.27267
total cost =============================: 4.95969
total  =================================: 186
8.24274
```



下面开始把string换掉了：

整个项目中的字符串操作大体上就是先解析读进来的char*，找到start和end之后构造s=string，把s传进具体的处理函数使用，过程中出现了对s的slect和cut，slect即筛选，有的质量很低的s就忽略掉了，cut即s.substr等，然后把处理之后的s拼接成一个大的ans，输出到文件中。

根据对string的具体了解，它的构造函数（以及substr）(start,end) or (start,len)都是把[s,e]的数据复制了一遍，而s.c_str()则是简单的指针指过来，ans+=s的操作

这个过程中可以优化的点就是，整个string其实就是一个中间变量，我们可以不要他，直接用start和end来表示s。但是直接这样大改要修改的地方太多了， 我们先从最后一个步骤优化起，即最后不拼接了。



### 0114

好家伙，咕了一个月了，回家就是吃吃吃睡睡睡，摸了4天鱼了，不能再摸了，开始干活了。

关于mDumplite的bug还是没有解决，先去换string吧。

./rabbit_qc -w 1 -i ../data/test.fq -o p.fq --umi --umi_loc index1

开启umi，但是依旧不是热点，而且cut里面的substr似乎并不慢（初步猜测是s=s.substr的原因）。

### 0115

./rabbit_qc -w 1 -i ../data/test.fq -o p.fq --umi --umi_loc index1 index2 -t 4 -5 -3 -g -x

开启PolyX，也不算是热点。

还差一个AdapterTrimmer::trimBySequence没有打开，带adapter的数据正在下，下午开开试试。

目前se的所有功能似乎都没有删除的操作，对字符串的修改都是cut，因此string换char*的思路是可以的。

研究一下项目的输出系统：

首先，mOptions->split.enabled即，输出要不要split成多个文件，默认是关闭的；

initOutput()中声明了mLeftWriter，默认情况下就是new了1e7的数组，开了一个输出流；

每一个Task（即处理一个Pack）之后都会把拼接好的临时变量outstr复制到一个新的空间，然后把地址和大小存好，准备最后统一输出；

![IMG_0480.PNG](https://i.loli.net/2021/01/15/CIE6swqjOWuGDYU.png)



中间data的变化大体就是这样的，初步的思路是Read中不再有新的string，而是两个指针。

晚上找zz讨论了一下，他说这样有点激进，那就先把outstr去掉吧，还能剩下Read类：

|            | test.fq | SRR2496709_1.fastq |
| ---------- | ------- | ------------------ |
| STD        | 7.30    | 48.8               |
| Now        | 6.31    | 43.0               |
| no out_str | 5.94    | 38.58              |

但是这个版本没有及时释放空间，内存需要从1g变成了6g(3.5g输入)，6g显然不合适，只能加上delete。

关于空间上的问题，似乎还没有那么简单：


- [x] 每次循环的r1和or1是两个新的空间吗，还是两个指针指向一个空间（或者开不开-t结果不一样）
- [ ] Read的delete调用的是啥，会释放，如果会，那if (r1 != or1 && r1 != NULL)delete r1;会发生什么
- [ ] std::vector<Read *>().swap(pack->data);是个啥玩意，是因为之前data都delete过了，所以弄成空的，防止delete pack报错吗

经过测试，r1有时候可能为NULL，除此之外，r1和or1都是一样的，即不同的指针，指向同一块空间，因此后面的delete就需要改改；

关于delete没有具体的测试，但是基本上可以确定就是会释放，那个if 不会执行；

是的；

|            | test.fq       | SRR2496709_1.fastq | SRR2530740.sra.fastq |
| ---------- | ------------- | ------------------ | -------------------- |
| STD        | 5.50+1.65=8.7 | 40.46+10.18=56.81  | 68.60+20.58=92.88    |
| Now        | Xxx           | Xxx                | Xxx                  |
| no out_str | 3.56+1.62=7.0 | 32.26+8.92=46.88   | 49.54+20.68=79.38    |



```c++
./rabbit_qc -w 1 -i ../data/test.fq -o p.fq -U --umi_loc=read1 --umi_len=8
total getPreStats1()->statRead(or1) ====: 1.45745
total mDuplicate->statRead(or1) ========: 0.50807
total mOptions->indexFilter()  =========: 0.05273
total mUmiProcessor->process(or1) ======: 0.98117
total mFilter->trimAndCut() ============: 0.06615
total PolyX::trimPolyG() ===============: 0.98117
total trimBySequence ===================: 0.05326
total r1->resize() =====================: 0.05256
total mFilter->passFilter(r1) ==========: 0.14185
total addFilterResult(result) ==========: 0.05432
total outstr += r1->toString() =========: 0.23384
total getPostStats1()->statRead(r1) ====: 1.01119
total delete r1 ========================: 0.05133
total ready output ========================: 1.70261
total costTotel ========================: 4.71656
total cost =============================: 5.45134
total  =================================: 186
total format =================================: 2.82441
Read1 before filtering:
total reads: 2000000
total bases: 200000000
Q20 bases: 196873669(98.4368%)
Q30 bases: 192134023(96.067%)

Read1 after filtering:
total reads: 1977633
total bases: 181942236
Q20 bases: 180366667(99.134%)
Q30 bases: 175978730(96.7223%)

Filtering result:
reads passed filter: 1977633
reads failed due to low quality: 22143
reads failed due to too many N: 224
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

Duplication rate (may be overestimated since this is SE data): 0.153535%

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../data/test.fq -o p.fq -U --umi_loc=read1 --umi_len=8 
rabbit_qc v0.0.1, time used: 11.7879 seconds
```





#### 0327

gkd gkd

现在的版本：修改了mduplicate的bug，写线程还在原位置，加了check，如果没有指定-o，则不生成outstr等，指定了的话直接存指针，减少了一步内存拷贝。

#### 0329

on Mac thread 1:

```
➜  RabbitQCPlus git:(ylf) ✗ ./rabbit_qc -w 1 -i ../data/test.fq
8 CPUs detected
Detecting adapter sequence for read1...
No adapter detected for read1

mKeyLenInBase 12
producer.join
threads.join
total getPreStats1()->statRead(or1) ====: 1.03896
total mDuplicate->statRead(or1) ========: 0.61504
total mOptions->indexFilter()  =========: 0.06373
total mUmiProcessor->process(or1) ======: 0.06161
total mFilter->trimAndCut() ============: 0.06570
total PolyX::trimPolyG() ===============: 0.06161
total trimBySequence ===================: 0.06265
total r1->resize() =====================: 0.06206
total mFilter->passFilter(r1) ==========: 0.10780
total addFilterResult(result) ==========: 0.06479
total outstr += r1->toString() =========: 0.06240
total getPostStats1()->statRead(r1) ====: 0.98106
total delete r1 ========================: 0.93583
total ready output ========================: 0.00093
total costTotel ========================: 4.18400
total cost =============================: 5.02375
total  =================================: 186
total format =================================: 3.05091
Read1 before filtering:
total reads: 2000000
total bases: 200000000
Q20 bases: 196873669(98.4368%)
Q30 bases: 192134023(96.067%)

Read1 after filtering:
total reads: 1979553
total bases: 197955300
Q20 bases: 196261878(99.1445%)
Q30 bases: 191688154(96.8341%)

Filtering result:
reads passed filter: 1979553
reads failed due to low quality: 20180
reads failed due to too many N: 267
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

Duplication rate (may be overestimated since this is SE data): 0.153535%

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../data/test.fq
rabbit_qc v0.0.1, time used: 9.28793 seconds


➜  STD git:(STD) ✗ ./rabbit_qc -w 1 -i ../data/test.fq
8 CPUs detected
Detecting adapter sequence for read1...
No adapter detected for read1

all thread created
producer join
consumers join
writer join
total getPreStats1()->statRead(or1) ====: 1.08322
total mDuplicate->statRead(or1) ========: 1.70167
total mOptions->indexFilter()  =========: 0.06277
total mUmiProcessor->process(or1) ======: 0.06253
total mFilter->trimAndCut() ============: 0.06567
total PolyX::trimPolyG() ===============: 0.06253
total trimBySequence ===================: 0.06258
total r1->resize() =====================: 0.06189
total mFilter->passFilter(r1) ==========: 0.13293
total addFilterResult(result) ==========: 0.06650
total outstr += r1->toString() =========: 1.22935
total getPostStats1()->statRead(r1) ====: 0.91797
total delete r1 ========================: 0.96272
total ready output ========================: 0.00084
total costTotel ========================: 6.47184
total cost =============================: 7.34789
total  =================================: 186
total format =================================: 3.14067
Read1 before filtering:
total reads: 2000000
total bases: 200000000
Q20 bases: 196873669(98.4368%)
Q30 bases: 192134023(96.067%)

Read1 after filtering:
total reads: 1979553
total bases: 197955300
Q20 bases: 196261878(99.1445%)
Q30 bases: 191688154(96.8341%)

Filtering result:
reads passed filter: 1979553
reads failed due to low quality: 20180
reads failed due to too many N: 267
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

Duplication rate (may be overestimated since this is SE data): 0.153535%

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../data/test.fq
rabbit_qc v0.0.1, time used: 11.7177 seconds




➜  RabbitQCPlus git:(ylf) ✗ ./rabbit_qc -w 1 -i ../data/SRR2496709_1.fastq
8 CPUs detected
Detecting adapter sequence for read1...
Illumina TruSeq Adapter Read 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

mKeyLenInBase 12
producer.join
threads.join
total getPreStats1()->statRead(or1) ====: 9.58512
total mDuplicate->statRead(or1) ========: 1.89019
total mOptions->indexFilter()  =========: 0.40946
total mUmiProcessor->process(or1) ======: 0.39962
total mFilter->trimAndCut() ============: 0.42816
total PolyX::trimPolyG() ===============: 0.39962
total trimBySequence ===================: 9.00838
total r1->resize() =====================: 0.40214
total mFilter->passFilter(r1) ==========: 0.71852
total addFilterResult(result) ==========: 0.41392
total outstr += r1->toString() =========: 0.40031
total getPostStats1()->statRead(r1) ====: 6.21279
total delete r1 ========================: 6.47128
total ready output ========================: 0.00556
total costTotel ========================: 36.74017
total cost =============================: 42.22726
total  =================================: 1124
total format =================================: 18.93405
Read1 before filtering:
total reads: 12607412
total bases: 1260741163
Q20 bases: 1160816635(92.0741%)
Q30 bases: 548113298(43.4755%)

Read1 after filtering:
total reads: 12375424
total bases: 1229937278
Q20 bases: 1146407885(93.2086%)
Q30 bases: 543968479(44.2273%)

Filtering result:
reads passed filter: 12375424
reads failed due to low quality: 228241
reads failed due to too many N: 3687
reads failed due to too short: 60
reads with adapter trimmed: 481110
bases trimmed due to adapters: 7646566

Duplication rate (may be overestimated since this is SE data): 64.1312%

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../data/SRR2496709_1.fastq
rabbit_qc v0.0.1, time used: 62.2956 seconds


➜  STD git:(STD) ✗ ./rabbit_qc -w 1 -i ../data/SRR2496709_1.fastq
8 CPUs detected
Detecting adapter sequence for read1...
Illumina TruSeq Adapter Read 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

all thread created
producer join
consumers join
writer join
total getPreStats1()->statRead(or1) ====: 10.66436
total mDuplicate->statRead(or1) ========: 6.25781
total mOptions->indexFilter()  =========: 0.42511
total mUmiProcessor->process(or1) ======: 0.42186
total mFilter->trimAndCut() ============: 0.44781
total PolyX::trimPolyG() ===============: 0.42186
total trimBySequence ===================: 10.28335
total r1->resize() =====================: 0.41633
total mFilter->passFilter(r1) ==========: 0.92623
total addFilterResult(result) ==========: 0.44656
total outstr += r1->toString() =========: 9.26387
total getPostStats1()->statRead(r1) ====: 7.23739
total delete r1 ========================: 6.78670
total ready output ========================: 0.00575
total costTotel ========================: 53.99258
total cost =============================: 59.95749
total  =================================: 1124
total format =================================: 21.38945
Read1 before filtering:
total reads: 12607412
total bases: 1260741163
Q20 bases: 1160816627(92.0741%)
Q30 bases: 548113277(43.4755%)

Read1 after filtering:
total reads: 12375423
total bases: 1229937215
Q20 bases: 1146407842(93.2086%)
Q30 bases: 543968436(44.2273%)

Filtering result:
reads passed filter: 12375423
reads failed due to low quality: 228242
reads failed due to too many N: 3687
reads failed due to too short: 60
reads with adapter trimmed: 481110
bases trimmed due to adapters: 7646566

Duplication rate (may be overestimated since this is SE data): 64.1312%

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../data/SRR2496709_1.fastq
rabbit_qc v0.0.1, time used: 82.6602 seconds
```



on 6148 thread 1:

```
ylf@gold6148:~/QC/STD$ ./rm.sh && time ./rabbit_qc -w 1 -i ../../data/test.fq
40 CPUs detected
Detecting adapter sequence for read1...
No adapter detected for read1

all thread created
producer join
consumers join
writer join
total getPreStats1()->statRead(or1) ====: 0.98390
total mDuplicate->statRead(or1) ========: 1.37794
total mOptions->indexFilter()  =========: 0.04385
total mUmiProcessor->process(or1) ======: 0.10784
total mFilter->trimAndCut() ============: 0.05045
total PolyX::trimPolyG() ===============: 0.10784
total trimBySequence ===================: 0.04439
total r1->resize() =====================: 0.04411
total mFilter->passFilter(r1) ==========: 0.11005
total addFilterResult(result) ==========: 0.04653
total outstr += r1->toString() =========: 0.46217
total getPostStats1()->statRead(r1) ====: 0.93202
total delete r1 ========================: 0.23774
total ready output ========================: 0.01627
total costTotel ========================: 4.48514
total cost =============================: 5.10614
total  =================================: 186
total format =================================: 1.43115
Read1 before filtering:
total reads: 2000000
total bases: 200000000
Q20 bases: 196873669(98.4368%)
Q30 bases: 192134023(96.067%)

Read1 after filtering:
total reads: 1979553
total bases: 197955300
Q20 bases: 196261878(99.1445%)
Q30 bases: 191688154(96.8341%)

Filtering result:
reads passed filter: 1979553
reads failed due to low quality: 20180
reads failed due to too many N: 267
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

Duplication rate (may be overestimated since this is SE data): 0.153535%

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../../data/test.fq
rabbit_qc v0.0.1, time used: 7.14761 seconds

real	0m7.156s
user	0m9.127s
sys	0m0.482s

ylf@gold6148:~/QC/RabbitQC$ ./rm.sh && time ./rabbit_qc -w 1 -i ../../data/test.fq
40 CPUs detected
Detecting adapter sequence for read1...
No adapter detected for read1

mKeyLenInBase 12
producer.join
threads.join
total getPreStats1()->statRead(or1) ====: 0.94659
total mDuplicate->statRead(or1) ========: 0.68595
total mOptions->indexFilter()  =========: 0.04383
total mUmiProcessor->process(or1) ======: 0.04287
total mFilter->trimAndCut() ============: 0.04946
total PolyX::trimPolyG() ===============: 0.04287
total trimBySequence ===================: 0.04391
total r1->resize() =====================: 0.04414
total mFilter->passFilter(r1) ==========: 0.09220
total addFilterResult(result) ==========: 0.04671
total outstr += r1->toString() =========: 0.04407
total getPostStats1()->statRead(r1) ====: 0.87878
total delete r1 ========================: 0.19308
total ready output ========================: 0.08632
total costTotel ========================: 3.15576
total cost =============================: 3.75207
total  =================================: 186
total format =================================: 1.46433
Read1 before filtering:
total reads: 2000000
total bases: 200000000
Q20 bases: 196873669(98.4368%)
Q30 bases: 192134023(96.067%)

Read1 after filtering:
total reads: 1979553
total bases: 197955300
Q20 bases: 196261878(99.1445%)
Q30 bases: 191688154(96.8341%)

Filtering result:
reads passed filter: 1979553
reads failed due to low quality: 20180
reads failed due to too many N: 267
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

Duplication rate (may be overestimated since this is SE data): 0.153535%

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../../data/test.fq
rabbit_qc v0.0.1, time used: 5.97684 seconds

real	0m5.986s
user	0m7.657s
sys	0m0.716s



ylf@gold6148:~/QC/STD$ ./rm.sh && time ./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq
40 CPUs detected
Detecting adapter sequence for read1...
No adapter detected for read1

all thread created
producer join
consumers join
writer join
total getPreStats1()->statRead(or1) ====: 13.42695
total mDuplicate->statRead(or1) ========: 19.12346
total mOptions->indexFilter()  =========: 0.60050
total mUmiProcessor->process(or1) ======: 1.07507
total mFilter->trimAndCut() ============: 0.71319
total PolyX::trimPolyG() ===============: 1.07507
total trimBySequence ===================: 0.60844
total r1->resize() =====================: 0.60799
total mFilter->passFilter(r1) ==========: 1.54834
total addFilterResult(result) ==========: 0.63633
total outstr += r1->toString() =========: 6.30583
total getPostStats1()->statRead(r1) ====: 12.78121
total delete r1 ========================: 3.28333
total ready output ========================: 0.21955
total costTotel ========================: 61.32007
total cost =============================: 69.85758
total  =================================: 2552
total format =================================: 18.90963
Read1 before filtering:
total reads: 27497479
total bases: 2749747900
Q20 bases: 2707938680(98.4795%)
Q30 bases: 2645554729(96.2108%)

Read1 after filtering:
total reads: 27219050
total bases: 2721905000
Q20 bases: 2699552056(99.1788%)
Q30 bases: 2639324248(96.9661%)

Filtering result:
reads passed filter: 27219050
reads failed due to low quality: 274717
reads failed due to too many N: 3712
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

Duplication rate (may be overestimated since this is SE data): 0.159583%

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq
rabbit_qc v0.0.1, time used: 89.6717 seconds

real	1m29.680s
user	1m31.963s
sys	0m3.520s


ylf@gold6148:~/QC/RabbitQC$ ./rm.sh && time ./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq 40 CPUs detected
Detecting adapter sequence for read1...
No adapter detected for read1

mKeyLenInBase 12
producer.join
threads.join
total getPreStats1()->statRead(or1) ====: 12.76642
total mDuplicate->statRead(or1) ========: 6.90075
total mOptions->indexFilter()  =========: 0.59985
total mUmiProcessor->process(or1) ======: 0.59829
total mFilter->trimAndCut() ============: 0.74826
total PolyX::trimPolyG() ===============: 0.59829
total trimBySequence ===================: 0.60584
total r1->resize() =====================: 0.60534
total mFilter->passFilter(r1) ==========: 1.26489
total addFilterResult(result) ==========: 0.61683
total outstr += r1->toString() =========: 0.60175
total getPostStats1()->statRead(r1) ====: 12.22309
total delete r1 ========================: 2.64687
total ready output ========================: 1.03176
total costTotel ========================: 40.78560
total cost =============================: 48.94914
total  =================================: 2552
total format =================================: 17.36762

Read1 before filtering:
total reads: 27497479
total bases: 2749747900
Q20 bases: 2707938680(98.4795%)
Q30 bases: 2645554729(96.2108%)

Read1 after filtering:
total reads: 27219050
total bases: 2721905000
Q20 bases: 2699552056(99.1788%)
Q30 bases: 2639324248(96.9661%)

Filtering result:
reads passed filter: 27219050
reads failed due to low quality: 274717
reads failed due to too many N: 3712
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

Duplication rate (may be overestimated since this is SE data): 0.159583%

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq
rabbit_qc v0.0.1, time used: 67.9869 seconds

real	1m7.997s
user	1m10.412s
sys	0m3.542s



ylf@gold6148:~/QC/STD$ ./rm.sh && time ./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq -o p.fq
40 CPUs detected
Detecting adapter sequence for read1...
No adapter detected for read1

all thread created
producer join
consumers join
writer join
total getPreStats1()->statRead(or1) ====: 13.40239
total mDuplicate->statRead(or1) ========: 19.18981
total mOptions->indexFilter()  =========: 0.60235
total mUmiProcessor->process(or1) ======: 1.09326
total mFilter->trimAndCut() ============: 0.67057
total PolyX::trimPolyG() ===============: 1.09326
total trimBySequence ===================: 0.60988
total r1->resize() =====================: 0.60936
total mFilter->passFilter(r1) ==========: 1.52353
total addFilterResult(result) ==========: 0.64401
total outstr += r1->toString() =========: 6.78953
total getPostStats1()->statRead(r1) ====: 13.11158
total delete r1 ========================: 3.33387
total ready output ========================: 1.02309
total costTotel ========================: 62.18672
total cost =============================: 70.72930
total  =================================: 2552
total format =================================: 19.04990
Read1 before filtering:
total reads: 27497479
total bases: 2749747900
Q20 bases: 2707938680(98.4795%)
Q30 bases: 2645554729(96.2108%)

Read1 after filtering:
total reads: 27219050
total bases: 2721905000
Q20 bases: 2699552056(99.1788%)
Q30 bases: 2639324248(96.9661%)

Filtering result:
reads passed filter: 27219050
reads failed due to low quality: 274717
reads failed due to too many N: 3712
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

Duplication rate (may be overestimated since this is SE data): 0.159583%

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq -o p.fq
rabbit_qc v0.0.1, time used: 91.4922 seconds

real	1m31.504s
user	1m34.372s
sys	0m8.842s


ylf@gold6148:~/QC/RabbitQC$ ./rm.sh && time ./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq  -o p.fq
40 CPUs detected
Detecting adapter sequence for read1...
No adapter detected for read1

mKeyLenInBase 12
producer.join
threads.join
leftWriterThread->join
total getPreStats1()->statRead(or1) ====: 13.27610
total mDuplicate->statRead(or1) ========: 5.24192
total mOptions->indexFilter()  =========: 0.60064
total mUmiProcessor->process(or1) ======: 1.14907
total mFilter->trimAndCut() ============: 0.72638
total PolyX::trimPolyG() ===============: 1.14907
total trimBySequence ===================: 0.60698
total r1->resize() =====================: 0.60780
total mFilter->passFilter(r1) ==========: 1.30622
total addFilterResult(result) ==========: 0.62687
total outstr += r1->toString() =========: 0.64608
total getPostStats1()->statRead(r1) ====: 12.08483
total delete r1 ========================: 0.60040
total ready output ========================: 6.29991
total costTotel ========================: 38.08050
total cost =============================: 46.26352
total  =================================: 2552
total format =================================: 19.72727
Read1 before filtering:
total reads: 27497479
total bases: 2749747900
Q20 bases: 2707938680(98.4795%)
Q30 bases: 2645554729(96.2108%)

Read1 after filtering:
total reads: 27219050
total bases: 2721905000
Q20 bases: 2699552056(99.1788%)
Q30 bases: 2639324248(96.9661%)

Filtering result:
reads passed filter: 27219050
reads failed due to low quality: 274717
reads failed due to too many N: 3712
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

Duplication rate (may be overestimated since this is SE data): 0.159583%

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq -o p.fq
rabbit_qc v0.0.1, time used: 72.9528 seconds

real	1m12.968s
user	1m15.486s
sys	0m8.658s
```

下面研究优化state：

首先比较容易想到的是对统计信息的for循环做向量化，但是由于访存的不规律，不太好弄。暂时先把原来的[base] [pos]换成[pos] [base]提高cache的命中率，并且方便后面的向量化：

emmmm没啥用啊，如果把12维调换，然后8个8个的来处理，看汇编发现并没有自动向量化，毕竟内存中本来就不连续，-Ofast什么的只是把循环展开了，关于cache命中率可能也没啥用。

现在的版本把12维调换了（暂时没什么用，可以先保留着），然后把一维的两个变量摘出来做自动向量化。

```
ylf@gold6148:~/QC/RabbitQC$ ./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq
40 CPUs detected
Detecting adapter sequence for read1...
No adapter detected for read1

mKeyLenInBase 12
producer.join
threads.join
total getPreStats1()->statRead(or1) ====: 12.34479
total mDuplicate->statRead(or1) ========: 4.93429
total mOptions->indexFilter()  =========: 0.59897
total mUmiProcessor->process(or1) ======: 0.89748
total mFilter->trimAndCut() ============: 0.74583
total PolyX::trimPolyG() ===============: 0.89748
total trimBySequence ===================: 0.60313
total r1->resize() =====================: 0.60471
total mFilter->passFilter(r1) ==========: 1.27222
total addFilterResult(result) ==========: 0.60432
total outstr += r1->toString() =========: 0.59924
total getPostStats1()->statRead(r1) ====: 11.93036
total delete r1 ========================: 2.69691
total ready output ========================: 1.12396
total costTotel ========================: 38.43762
total cost =============================: 46.54946
total  =================================: 2552
total format =================================: 17.74351
Read1 before filtering:
total reads: 27497479
total bases: 2749747900
Q20 bases: 2707938680(98.4795%)
Q30 bases: 2645554729(96.2108%)

Read1 after filtering:
total reads: 27219050
total bases: 2721905000
Q20 bases: 2699552056(99.1788%)
Q30 bases: 2639324248(96.9661%)

Filtering result:
reads passed filter: 27219050
reads failed due to low quality: 274717
reads failed due to too many N: 3712
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

Duplication rate (may be overestimated since this is SE data): 0.159583%

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq
rabbit_qc v0.0.1, time used: 66.4818 seconds
```

#### 0330

state的向量化暂时没啥思路了，唯一的可能性就是*8/16。

format中new read时有一遍多余的内存拷贝，去掉大约能快5s(SRR2530740.sra.fastq 7.5G)。

#### 0331

手动向量化，暂时处理A[100]的数据，测试了set和load的速度区别

```
_mm512_set_epi64 : 66.496
_mm512_load_epi64 :  65.8957
```

#### 0402

现在有点子奇怪，duplicate的统计函数，-march=native能比-mavx2满好几秒，现在看看这两段的汇编找找原因。

好像是不知道为啥native的版本没有向量化，手写了一部分，现在总时间和全部自动向量化差不多。。。。

....

添加了mm256的版本，但是还有点bug。

#### 0407

现在自动向量化的版本

```
ylf@gold6148:~/QC/RabbitQC$ ./rm.sh && time ./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq
40 CPUs detected
Detecting adapter sequence for read1...
No adapter detected for read1

mKeyLenInBase 12
producer.join
threads.join
total getPreStats1()->statRead(or1) ====: 12.36815
total mDuplicate->statRead(or1) ========: 4.96933
total mOptions->indexFilter()  =========: 0.59926
total mUmiProcessor->process(or1) ======: 0.95765
total mFilter->trimAndCut() ============: 0.66535
total PolyX::trimPolyG() ===============: 0.95765
total trimBySequence ===================: 0.60368
total r1->resize() =====================: 0.60566
total mFilter->passFilter(r1) ==========: 1.24261
total addFilterResult(result) ==========: 0.65178
total outstr += r1->toString() =========: 0.60243
total getPostStats1()->statRead(r1) ====: 11.97133
total delete r1 ========================: 2.70405
total ready output ========================: 1.22500
total costTotel ========================: 38.54645
total cost =============================: 46.69251
total  =================================: 2552
total format =================================: 13.03277
Read1 before filtering:
total reads: 27497479
total bases: 2749747900
Q20 bases: 2707938680(98.4795%)
Q30 bases: 2645554729(96.2108%)

Read1 after filtering:
total reads: 27219050
total bases: 2721905000
Q20 bases: 2699552056(99.1788%)
Q30 bases: 2639324248(96.9661%)

Filtering result:
reads passed filter: 27219050
reads failed due to low quality: 274717
reads failed due to too many N: 3712
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

Duplication rate (may be overestimated since this is SE data): 0.159583%

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq
rabbit_qc v0.0.1, time used: 61.641 seconds

real	1m1.681s
user	1m3.682s
sys	0m4.101s
```

把state里面的long换成unsigned int。统计信息。

```
ylf@gold6148:~/QC/STD$ ./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq -o p.fq
40 CPUs detected
Detecting adapter sequence for read1...
No adapter detected for read1

all thread created
producer join
consumers join
writer join
total getPreStats1()->statRead(or1) ====: 13.09763
total mDuplicate->statRead(or1) ========: 19.80397
total mOptions->indexFilter()  =========: 0.60315
total mUmiProcessor->process(or1) ======: 1.08190
total mFilter->trimAndCut() ============: 0.68614
total PolyX::trimPolyG() ===============: 1.08190
total trimBySequence ===================: 0.60624
total r1->resize() =====================: 0.60851
total mFilter->passFilter(r1) ==========: 1.52112
total addFilterResult(result) ==========: 0.66938
total outstr += r1->toString() =========: 6.94076
total getPostStats1()->statRead(r1) ====: 13.00677
total delete r1 ========================: 3.29540
total ready output ========================: 0.84695
total costTotel ========================: 62.52968
total cost =============================: 71.07536
total  =================================: 2552
total format =================================: 18.83949
Read1 before filtering:
total reads: 27497479
total bases: 2749747900
Q20 bases: 2707938680(98.4795%)
Q30 bases: 2645554729(96.2108%)

Read1 after filtering:
total reads: 27219050
total bases: 2721905000
Q20 bases: 2699552056(99.1788%)
Q30 bases: 2639324248(96.9661%)

Filtering result:
reads passed filter: 27219050
reads failed due to low quality: 274717
reads failed due to too many N: 3712
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

Duplication rate (may be overestimated since this is SE data): 0.159583%

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq -o p.fq
rabbit_qc v0.0.1, time used: 91.6184 seconds


ylf@gold6148:~/QC/RabbitQC$ ./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq -o p.fq
40 CPUs detected
Detecting adapter sequence for read1...
No adapter detected for read1

mKeyLenInBase 12
producer.join
threads.join
leftWriterThread->join
total getPreStats1()->statRead(or1) ====: 10.84282
total mDuplicate->statRead(or1) ========: 5.22429
total mOptions->indexFilter()  =========: 0.60391
total mUmiProcessor->process(or1) ======: 1.13022
total mFilter->trimAndCut() ============: 0.67868
total PolyX::trimPolyG() ===============: 1.13022
total trimBySequence ===================: 0.60568
total r1->resize() =====================: 0.60597
total mFilter->passFilter(r1) ==========: 1.24626
total addFilterResult(result) ==========: 0.61152
total outstr += r1->toString() =========: 0.65207
total getPostStats1()->statRead(r1) ====: 9.93785
total delete r1 ========================: 0.60108
total ready output ========================: 6.51360
total costTotel ========================: 33.34706
total cost =============================: 41.48986
total  =================================: 2552
total format =================================: 14.99609

Read1 before filtering:
total reads: 27497479
total bases: 2749747900
Q20 bases: 2707938680(98.4795%)
Q30 bases: 2645554729(96.2108%)

Read1 after filtering:
total reads: 27219050
total bases: 2721905000
Q20 bases: 2699552056(99.1788%)
Q30 bases: 2639324248(96.9661%)

Filtering result:
reads passed filter: 27219050
reads failed due to low quality: 274717
reads failed due to too many N: 3712
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

Duplication rate (may be overestimated since this is SE data): 0.159583%

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq -o p.fq
rabbit_qc v0.0.1, time used: 63.8425 seconds






ylf@gold6148:~/QC/STD$ ./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq
40 CPUs detected
Detecting adapter sequence for read1...
No adapter detected for read1

all thread created
producer join
consumers join
writer join
total getPreStats1()->statRead(or1) ====: 13.29855
total mDuplicate->statRead(or1) ========: 20.11957
total mOptions->indexFilter()  =========: 0.60335
total mUmiProcessor->process(or1) ======: 1.12197
total mFilter->trimAndCut() ============: 0.66747
total PolyX::trimPolyG() ===============: 1.12197
total trimBySequence ===================: 0.60532
total r1->resize() =====================: 0.60786
total mFilter->passFilter(r1) ==========: 1.52133
total addFilterResult(result) ==========: 0.66952
total outstr += r1->toString() =========: 6.27320
total getPostStats1()->statRead(r1) ====: 12.67426
total delete r1 ========================: 3.29680
total ready output ========================: 0.23269
total costTotel ========================: 62.06691
total cost =============================: 70.61276
total  =================================: 2552
total format =================================: 19.12414
Read1 before filtering:
total reads: 27497479
total bases: 2749747900
Q20 bases: 2707938680(98.4795%)
Q30 bases: 2645554729(96.2108%)

Read1 after filtering:
total reads: 27219050
total bases: 2721905000
Q20 bases: 2699552056(99.1788%)
Q30 bases: 2639324248(96.9661%)

Filtering result:
reads passed filter: 27219050
reads failed due to low quality: 274717
reads failed due to too many N: 3712
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

Duplication rate (may be overestimated since this is SE data): 0.159583%

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq
rabbit_qc v0.0.1, time used: 90.7219 seconds

ylf@gold6148:~/QC/RabbitQC$ ./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq
40 CPUs detected
Detecting adapter sequence for read1...
No adapter detected for read1

mKeyLenInBase 12
producer.join
threads.join
total getPreStats1()->statRead(or1) ====: 10.22650
total mDuplicate->statRead(or1) ========: 5.10257
total mOptions->indexFilter()  =========: 0.60111
total mUmiProcessor->process(or1) ======: 1.11235
total mFilter->trimAndCut() ============: 0.68214
total PolyX::trimPolyG() ===============: 1.11235
total trimBySequence ===================: 0.60651
total r1->resize() =====================: 0.60762
total mFilter->passFilter(r1) ==========: 1.26068
total addFilterResult(result) ==========: 0.61204
total outstr += r1->toString() =========: 0.60256
total getPostStats1()->statRead(r1) ====: 9.83409
total delete r1 ========================: 2.67072
total ready output ========================: 1.22889
total costTotel ========================: 34.52398
total cost =============================: 42.65211
total  =================================: 2552
total format =================================: 12.99748
Read1 before filtering:
total reads: 27497479
total bases: 2749747900
Q20 bases: 2707938680(98.4795%)
Q30 bases: 2645554729(96.2108%)

Read1 after filtering:
total reads: 27219050
total bases: 2721905000
Q20 bases: 2699552056(99.1788%)
Q30 bases: 2639324248(96.9661%)

Filtering result:
reads passed filter: 27219050
reads failed due to low quality: 274717
reads failed due to too many N: 3712
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

Duplication rate (may be overestimated since this is SE data): 0.159583%

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq
rabbit_qc v0.0.1, time used: 57.6059 seconds
```

不打timer

```
ylf@gold6148:~/QC/STD$ ./run.sh
40 CPUs detected
Detecting adapter sequence for read1...
No adapter detected for read1

all thread created
producer join
consumers join
writer join
Read1 before filtering:
total reads: 27497479
total bases: 2749747900
Q20 bases: 2707938680(98.4795%)
Q30 bases: 2645554729(96.2108%)

Read1 after filtering:
total reads: 27219050
total bases: 2721905000
Q20 bases: 2699552056(99.1788%)
Q30 bases: 2639324248(96.9661%)

Filtering result:
reads passed filter: 27219050
reads failed due to low quality: 274717
reads failed due to too many N: 3712
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

Duplication rate (may be overestimated since this is SE data): 0.159583%

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq
rabbit_qc v0.0.1, time used: 75.0971 seconds


ylf@gold6148:~/QC/RabbitQC$ ./run.sh
40 CPUs detected
Detecting adapter sequence for read1...
No adapter detected for read1

mKeyLenInBase 12
producer.join
threads.join
Read1 before filtering:
total reads: 27497479
total bases: 2749747900
Q20 bases: 2707938680(98.4795%)
Q30 bases: 2645554729(96.2108%)

Read1 after filtering:
total reads: 27219050
total bases: 2721905000
Q20 bases: 2699552056(99.1788%)
Q30 bases: 2639324248(96.9661%)

Filtering result:
reads passed filter: 27219050
reads failed due to low quality: 274717
reads failed due to too many N: 3712
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

Duplication rate (may be overestimated since this is SE data): 0.159583%

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq
rabbit_qc v0.0.1, time used: 42.1073 seconds

real	0m42.120s
user	0m44.232s
sys	0m3.425s
```

开始弄弄pe的数据了

其实之前的优化基本上se pe的都通用，只有对于out_str的删除不一样，现在暂时先没有加这个，因为后期可能改的更激进一些，就不花时间在这里加了。下面记录了目前的运行时间，用的是1G的数据，但是两个大数据好像每次运行q20base都不一样。

```
ylf@gold6148:~/QC/STD$ ./rabbit_qc -w 1 -i ../../data/pe/1G_1.fastq -I ../../data/pe/1G_2.fastq
40 CPUs detected
parser cmd cost 0.000223875 seconds
adapter detected cost 0.0027709 seconds
total getPreStats1()->statRead(or1) ====: 5.24185
total mDuplicate->statRead(or1) ========: 0.74208
total mOptions->indexFilter()  =========: 0.08169
total mUmiProcessor->process(or1) ======: 0.09061
total mFilter->trimAndCut() ============: 0.10659
total PolyX::trimPolyG() ===============: 0.09061
total trimBySequence ===================: 6.98369
total r1->resize() =====================: 0.08317
total mFilter->passFilter(r1) ==========: 0.34177
total addFilterResult(result) ==========: 0.08671
total outstr += r1->toString() =========: 2.59094
total getPostStats1()->statRead(r1) ====: 3.86829
total delete r1 ========================: 0.72572
total ready output ========================: 0.22209
total costTotel ========================: 21.02912
total cost =============================: 22.20824
total  =================================: 257
total format =================================: 6.53803
。。。。
//uint
total getPreStats1()->statRead(or1) ====: 2.83674
//long
total getPreStats1()->statRead(or1) ====: 3.42496
//long vec512
total getPreStats1()->statRead(or1) ====: 3.48177
total mDuplicate->statRead(or1) ========: 0.29832
total mOptions->indexFilter()  =========: 0.08186
total mUmiProcessor->process(or1) ======: 0.09000
total mFilter->trimAndCut() ============: 0.11927
total PolyX::trimPolyG() ===============: 0.09000
total trimBySequence ===================: 5.57984
total r1->resize() =====================: 0.08256
total mFilter->passFilter(r1) ==========: 0.30770
total addFilterResult(result) ==========: 0.08412
total outstr += r1->toString() =========: 2.58373
total getPostStats1()->statRead(r1) ====: 2.56699
total delete r1 ========================: 0.71345
total ready output ========================: 0.22941
total costTotel ========================: 15.43081
total cost =============================: 16.56980
total  =================================: 257
total format =================================: 4.26026
。。。。
Read1 before filtering:
total reads: 3743702
total bases: 374370200
Q20 bases: 344779886(92.096%)
Q30 bases: 163671714(43.7192%)

Read1 after filtering:
total reads: 3491901
total bases: 346848212
Q20 bases: 327012029(94.281%)
Q30 bases: 158803010(45.7846%)

Read2 before filtering:
total reads: 3743702
total bases: 374370200
Q20 bases: 326390053(87.1838%)
Q30 bases: 127851348(34.151%)

Read2 aftering filtering:
total reads: 3491901
total bases: 346848212
Q20 bases: 315125942(90.8541%)
Q30 bases: 125820972(36.2755%)

Filtering result:
reads passed filter: 6983802
reads failed due to low quality: 501448
reads failed due to too many N: 2154
reads failed due to too short: 0
reads with adapter trimmed: 340622
bases trimmed due to adapters: 4733690

Duplication rate: 64.3233%

Insert size peak (evaluated by paired-end reads): 133

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../../data/pe/1G_1.fastq -I ../../data/pe/1G_2.fastq
rabbit_qc v0.0.1, time used: 29.244 seconds


ylf@gold6148:~/QC/RabbitQC$ ./rabbit_qc -w 1 -i ../../data/pe/1G_1.fastq -I ../../data/pe/1G_2.fastq
40 CPUs detected
parser cmd cost 0.000212908 seconds
adapter detected cost 0.00281096 seconds
mKeyLenInBase 12
finalPostStats1->mBufLen  100
total getPreStats1()->statRead(or1) ====: 2.83674
total mDuplicate->statRead(or1) ========: 0.29832
total mOptions->indexFilter()  =========: 0.08186
total mUmiProcessor->process(or1) ======: 0.09000
total mFilter->trimAndCut() ============: 0.11927
total PolyX::trimPolyG() ===============: 0.09000
total trimBySequence ===================: 5.57984
total r1->resize() =====================: 0.08256
total mFilter->passFilter(r1) ==========: 0.30770
total addFilterResult(result) ==========: 0.08412
total outstr += r1->toString() =========: 2.58373
total getPostStats1()->statRead(r1) ====: 2.56699
total delete r1 ========================: 0.71345
total ready output ========================: 0.22941
total costTotel ========================: 15.43081
total cost =============================: 16.56980
total  =================================: 257
total format =================================: 4.26026
Read1 before filtering:
total reads: 3743702
total bases: 374370200
Q20 bases: 344779886(92.096%)
Q30 bases: 163671714(43.7192%)

Read1 after filtering:
total reads: 3491901
total bases: 346848212
Q20 bases: 327012029(94.281%)
Q30 bases: 158803010(45.7846%)

Read2 before filtering:
total reads: 3743702
total bases: 374370200
Q20 bases: 326390053(87.1838%)
Q30 bases: 127851348(34.151%)

Read2 aftering filtering:
total reads: 3491901
total bases: 346848212
Q20 bases: 315125942(90.8541%)
Q30 bases: 125820972(36.2755%)

Filtering result:
reads passed filter: 6983802
reads failed due to low quality: 501448
reads failed due to too many N: 2154
reads failed due to too short: 0
reads with adapter trimmed: 340622
bases trimmed due to adapters: 4733690

Duplication rate: 64.3233%

Insert size peak (evaluated by paired-end reads): 133

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../../data/pe/1G_1.fastq -I ../../data/pe/1G_2.fastq
rabbit_qc v0.0.1, time used: 21.336 seconds



```

现在uint的情况下使用vec512进行测试

```
没啥用。。
```

把两个一维的统计信息直接留到所有的readPack处理完之后一块弄，快一丢丢

```
40 CPUs detected
parser cmd cost 0.000169992 seconds
Detecting adapter sequence for read1...
No adapter detected for read1

adapter detected cost 0.401944 seconds
mKeyLenInBase 12
producer.join
threads.join
total getPreStats1()->statRead(or1) ====: 9.71836
total mDuplicate->statRead(or1) ========: 4.97022
total mOptions->indexFilter()  =========: 0.60722
total mUmiProcessor->process(or1) ======: 1.08869
total mFilter->trimAndCut() ============: 0.68315
total PolyX::trimPolyG() ===============: 1.08869
total trimBySequence ===================: 0.60953
total r1->resize() =====================: 0.60990
total mFilter->passFilter(r1) ==========: 1.25422
total addFilterResult(result) ==========: 0.61375
total outstr += r1->toString() =========: 0.60632
total getPostStats1()->statRead(r1) ====: 9.28458
total delete r1 ========================: 2.65708
total ready output ========================: 1.22979
total costTotel ========================: 33.31543
total cost =============================: 41.51545
total  =================================: 2552
total format =================================: 13.00359
Read1 before filtering:
total reads: 27497479
total bases: 2749747900
Q20 bases: 2707938680(98.4795%)
Q30 bases: 2645554729(96.2108%)

Read1 after filtering:
total reads: 27219050
total bases: 2721905000
Q20 bases: 2699552056(99.1788%)
Q30 bases: 2639324248(96.9661%)

Filtering result:
reads passed filter: 27219050
reads failed due to low quality: 274717
reads failed due to too many N: 3712
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

Duplication rate (may be overestimated since this is SE data): 0.159583%

JSON report: RabbitQC.json
HTML report: RabbitQC.html

./rabbit_qc -w 1 -i ../../data/SRR2530740.sra.fastq
rabbit_qc v0.0.1, time used: 56.3852 seconds
```

<img src="/Users/ylf9811/Library/Containers/com.tencent.qq/Data/Library/Caches/Images/1F1A1E6E7FD34DD6E6E5553887E95D8F.jpg" alt="1F1A1E6E7FD34DD6E6E5553887E95D8F" style="zoom:50%;" />

<img src="/Users/ylf9811/Library/Containers/com.tencent.qq/Data/Library/Caches/Images/8AC851EB868DC59442BD46FD84A9CBCA.jpg" alt="8AC851EB868DC59442BD46FD84A9CBCA" style="zoom:50%;" />

就感觉自己最近挺搞笑的

修改bug之后，统计一下时间

```
init：rabbit_qc v0.0.1, time used: 75.0299 seconds
Vec512：rabbit_qc v0.0.1, time used: 40.981 seconds
Novec：rabbit_qc v0.0.1, time used: 41.2447 seconds
```



#### 0418

现在的程序测试一下时间：

首先，单端的数据 SRR2530740.sra.fastq（不带adapter） 不带输出  不包含adapter的检测

```
thread 1
STD ： 75.5552
Plus Vec512 + UseUint ： 40.7154
Plus Vec512 + UseLong ： 
Plus noVec  + UseLong ： 


thread 2
STD ： 41.0599
Plus Vec512 + UseUint ： 21.5165
Plus Vec512 + UseLong ： 
Plus noVec  + UseLong ： 


thread 4
STD ： 21.2547
Plus Vec512 + UseUint ： 11.0905
Plus Vec512 + UseLong ： 
Plus noVec  + UseLong ： 

thread 8
STD ： 11.5397
Plus Vec512 + UseUint ： 6.16717	
Plus Vec512 + UseLong ： 
Plus noVec  + UseLong ： 

thread 16
STD ： 6.1985
Plus Vec512 + UseUint ： 3.5977
Plus Vec512 + UseLong ： 
Plus noVec  + UseLong ： 
```

```
se 带adapter

thread 1
STD ： 46.307
Plus Vec512 + UseUint ： 22.7573
Plus Vec512 + UseLong ： 
Plus noVec  + UseLong ： 

thread 2
STD ： 23.2454
Plus Vec512 + UseUint ： 11.2436
Plus Vec512 + UseLong ： 
Plus noVec  + UseLong ： 

thread 4
STD ： 12.5085
Plus Vec512 + UseUint ： 5.95877
Plus Vec512 + UseLong ： 
Plus noVec  + UseLong ： 

thread 8
STD ： 6.62146
Plus Vec512 + UseUint ： 3.32433
Plus Vec512 + UseLong ： 
Plus noVec  + UseLong ： 

thread 16
STD ： 3.62546
Plus Vec512 + UseUint ： 1.91137
Plus Vec512 + UseLong ： 
Plus noVec  + UseLong ： 



pe 带adapter

thread 1
STD ： 93.2159
Plus Vec512 + UseUint ： 42.6254
Plus Vec512 + UseLong ： 
Plus noVec  + UseLong ： 

thread 2
STD ： 48.7778
Plus Vec512 + UseUint ： 21.0922
Plus Vec512 + UseLong ： 
Plus noVec  + UseLong ： 

thread 4
STD ： 24.3336
Plus Vec512 + UseUint ： 11.2543
Plus Vec512 + UseLong ： 
Plus noVec  + UseLong ： 

thread 8
STD ： 12.8506
Plus Vec512 + UseUint ： 6.13277
Plus Vec512 + UseLong ： 
Plus noVec  + UseLong ： 

thread 16
STD ： 6.87486
Plus Vec512 + UseUint ： 4.007
Plus Vec512 + UseLong ： 
Plus noVec  + UseLong ： 

```



