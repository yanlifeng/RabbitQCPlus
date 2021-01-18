## RabbitQC

TODO
  - [ ] 研究main函数中每一个参数的意义
  - [x] 打开trimAndCut
  - [x] 打开polg
  - [ ]  解决mDuplicate线程不安全问题
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
  - [ ] 为啥运行空间稳定1g不变，
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



中间data的变化答题就是这样的，初步的思路是Read中不再有新的string，而是两个指针。

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



