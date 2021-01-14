#include "seprocessor.h"
#include "fastqreader.h"
#include <iostream>
#include <unistd.h>
#include <functional>
#include <thread>
#include <memory.h>
#include "util.h"
#include "jsonreporter.h"
#include "htmlreporter.h"
#include "adaptertrimmer.h"
#include "polyx.h"


SingleEndProcessor::SingleEndProcessor(Options *opt) {
    mOptions = opt;
    mProduceFinished = false;
    mFinishedThreads = 0;
    mFilter = new Filter(opt);
//    printf("mFilter->test() start ...\n");
//    mFilter->test();
//    printf("mFilter->test() end ...");
    mOutStream = NULL;
    mZipFile = NULL;
    mUmiProcessor = new UmiProcessor(opt);
    mLeftWriter = NULL;

    mDuplicate = NULL;
    if (mOptions->duplicate.enabled) {
        mDuplicate = new Duplicate(mOptions);
    }

    //dsrc faster reader
    fastqPool = new dsrc::fq::FastqDataPool(128, 1 << 22);
}

SingleEndProcessor::~SingleEndProcessor() {
    delete mFilter;
    if (mDuplicate) {
        delete mDuplicate;
        mDuplicate = NULL;
    }
    //delete dsrc mem pool
    delete fastqPool;
}

void SingleEndProcessor::initOutput() {
    if (mOptions->out1.empty())
        return;
    mLeftWriter = new WriterThread(mOptions, mOptions->out1);
}

void SingleEndProcessor::closeOutput() {
    if (mLeftWriter) {
        delete mLeftWriter;
        mLeftWriter = NULL;
    }
}

void SingleEndProcessor::initConfig(ThreadConfig *config) {
    if (mOptions->out1.empty())
        return;

    if (mOptions->split.enabled) {
        config->initWriterForSplit();
    }
}

bool SingleEndProcessor::process() {
    if (!mOptions->split.enabled) {
        initOutput();
    }


    initPackRepository();
    std::thread
            producer(std::bind(&SingleEndProcessor::producerTask, this));

    //TODO: get the correct cycles
    int cycle = 151;
    ThreadConfig **configs = new ThreadConfig *[mOptions->thread];
    for (int t = 0; t < mOptions->thread; t++) {
        configs[t] = new ThreadConfig(mOptions, t, false);
        initConfig(configs[t]);
    }

    std::thread **threads = new thread *[mOptions->thread];
    for (int t = 0; t < mOptions->thread; t++) {
        threads[t] = new std::thread(std::bind(&SingleEndProcessor::consumerTask, this, configs[t]));
    }

    std::thread *leftWriterThread = NULL;
    if (mLeftWriter) {
//        printf("mLeftWriter is true, leftWriterThread pre ing ...\n");
        leftWriterThread = new std::thread(std::bind(&SingleEndProcessor::writeTask, this, mLeftWriter));
    }

    producer.join();
    for (int t = 0; t < mOptions->thread; t++) {
        threads[t]->join();
    }

    if (!mOptions->split.enabled) {
        if (leftWriterThread)
            leftWriterThread->join();
    }

    if (mOptions->verbose)
        loginfo("start to generate reports\n");

    // merge stats and read filter results
    vector<Stats *> preStats;
    vector<Stats *> postStats;
    vector<FilterResult *> filterResults;
    for (int t = 0; t < mOptions->thread; t++) {
        preStats.push_back(configs[t]->getPreStats1());
        postStats.push_back(configs[t]->getPostStats1());
        filterResults.push_back(configs[t]->getFilterResult());
    }
    Stats *finalPreStats = Stats::merge(preStats);
    Stats *finalPostStats = Stats::merge(postStats);
    FilterResult *finalFilterResult = FilterResult::merge(filterResults);


    // read filter results to the first thread's
    //TODO is this useful?

//    for (int t = 1; t < mOptions->thread; t++) {
//        preStats.push_back(configs[t]->getPreStats1());
//        postStats.push_back(configs[t]->getPostStats1());
//    }
    fstream out;

    fstream in1, in2, in3;
    string outFileName1 = "preStatsKmer";
    string outFileName2 = "postStatsKmer";
    string outFileName3 = "mDuplicateCount";
    string outFileName4 = "mDuplicateGC";
    string outFileName5 = "mDuplicateDups";

//
//    string checkGC = "../STD/mDuplicateGC";
//    string checkKmer = "../STD/mDuplicateDups";
//    string checkCnt = "../STD/mDuplicateCount";
//    // Read sampling file
//    uint8 *tempGC = new uint8[mDuplicate->mKeyLenInBit];
//    uint64 *tempKmer = new uint64[mDuplicate->mKeyLenInBit];
//    uint16 *tempCnt = new uint16[mDuplicate->mKeyLenInBit];
//    in1.open(checkGC, ios::in | ios::binary);
//    if (!in1) {
//        printf("Can't open file \"%s\"\n", checkGC.c_str());
//    } else {
//        in1.seekg(0, ios::beg);
//        in1.read(reinterpret_cast<char *>(tempGC), mDuplicate->mKeyLenInBit * sizeof(uint8));
//    }
//    in1.close();
//    in2.open(checkKmer, ios::in | ios::binary);
//    if (!in2) {
//        printf("Can't open file \"%s\"\n", checkKmer.c_str());
//    } else {
//        in2.seekg(0, ios::beg);
//        in2.read(reinterpret_cast<char *>(tempKmer), mDuplicate->mKeyLenInBit * sizeof(uint64));
//    }
//    in2.close();
//    in3.open(checkCnt, ios::in | ios::binary);
//    if (!in3) {
//        printf("Can't open file \"%s\"\n", checkCnt.c_str());
//    } else {
//        in3.seekg(0, ios::beg);
//        in3.read(reinterpret_cast<char *>(tempCnt), mDuplicate->mKeyLenInBit * sizeof(uint16));
//    }
//    in3.close();
//    printf("=================================================\n");
//    for (int i = 0; i < mDuplicate->mKeyLenInBit; i++) {
//        assert(tempCnt[i] == mDuplicate->mCounts[i]);
//        assert(tempKmer[i] == mDuplicate->mDups[i]);
//
////        for (int i = 0; i < 10000; i++) {
//        if (tempGC[i] != mDuplicate->mGC[i]) {
//            printf("GG on test %d  STD : %ld %d %d  Now : %ld %d %d\n", i, tempKmer[i], (int) tempGC[i], (int)
//                    tempCnt[i], mDuplicate->mDups[i], (int) mDuplicate->mGC[i], (int) mDuplicate->mCounts[i]);
//        }
////            printf("%ld\n", temp[i]);
//    }
//    printf("=================================================\n");

    out.open(outFileName1.c_str(), ios::out | ios::binary);
    out.seekp(0, ios::beg);
    out.write(reinterpret_cast<char *>(finalPreStats->mKmer), finalPreStats->mKmerBufLen * sizeof(long));
    out.close();

    out.open(outFileName2.c_str(), ios::out | ios::binary);
    out.seekp(0, ios::beg);
    out.write(reinterpret_cast<char *>(finalPostStats->mKmer), finalPostStats->mKmerBufLen * sizeof(long));
    out.close();

    out.open(outFileName3.c_str(), ios::out | ios::binary);
    out.seekp(0, ios::beg);
    out.write(reinterpret_cast<char *>(mDuplicate->mCounts), mDuplicate->mKeyLenInBit * sizeof(uint16));
    out.close();

    out.open(outFileName4.c_str(), ios::out | ios::binary);
    out.seekp(0, ios::beg);
    out.write(reinterpret_cast<char *>(mDuplicate->mGC), mDuplicate->mKeyLenInBit * sizeof(uint8));
    out.close();

    out.open(outFileName5.c_str(), ios::out | ios::binary);
    out.seekp(0, ios::beg);
    out.write(reinterpret_cast<char *>(mDuplicate->mDups), mDuplicate->mKeyLenInBit * sizeof(uint64));
    out.close();

    double cost = 0;
    double cost1 = 0;
    double cost2 = 0;
    double cost3 = 0;
    double cost4 = 0;
    double cost5 = 0;
    double cost6 = 0;
    double cost7 = 0;
    double cost8 = 0;
    double cost9 = 0;
    double cost10 = 0;
    double cost11 = 0;
    double cost12 = 0;
    double cost13 = 0;
    double costFormat = 0;
    int totCnt = 0;

    for (int t = 0; t < mOptions->thread; t++) {
        cost += configs[t]->cost;
        cost1 += configs[t]->cost1;
        cost2 += configs[t]->cost2;
        cost3 += configs[t]->cost3;
        cost4 += configs[t]->cost4;
        cost5 += configs[t]->cost5;
        cost6 += configs[t]->cost6;
        cost7 += configs[t]->cost7;
        cost8 += configs[t]->cost8;
        cost9 += configs[t]->cost9;
        cost10 += configs[t]->cost10;
        cost11 += configs[t]->cost11;
        cost12 += configs[t]->cost12;
        cost13 += configs[t]->cost13;
        totCnt += configs[t]->totCnt;
        costFormat += configs[t]->costFormat;
    }

    printf("total getPreStats1()->statRead(or1) ====: %.5f\n", cost1);
    printf("total mDuplicate->statRead(or1) ========: %.5f\n", cost2);
    printf("total mOptions->indexFilter()  =========: %.5f\n", cost3);
    printf("total mUmiProcessor->process(or1) ======: %.5f\n", cost4);
    printf("total mFilter->trimAndCut() ============: %.5f\n", cost5);
    printf("total PolyX::trimPolyG() ===============: %.5f\n", cost4);
    printf("total trimBySequence ===================: %.5f\n", cost7);
    printf("total r1->resize() =====================: %.5f\n", cost8);
    printf("total mFilter->passFilter(r1) ==========: %.5f\n", cost9);
    printf("total addFilterResult(result) ==========: %.5f\n", cost10);
    printf("total outstr += r1->toString() =========: %.5f\n", cost11);
    printf("total getPostStats1()->statRead(r1) ====: %.5f\n", cost12);
    printf("total delete r1 ========================: %.5f\n", cost13);
    printf("total costTotel ========================: %.5f\n",
           cost1 + cost2 + cost3 + cost4 + cost5 + cost6 + cost7 + cost8 + cost9 + cost10 + cost11 + cost12 + cost13);
    printf("total cost =============================: %.5f\n", cost);
    printf("total  =================================: %d\n", totCnt);
    printf("total format =================================: %.5f\n", costFormat);


    cerr << "Read1 before filtering:" << endl;
    finalPreStats->print();
    cerr << endl;
    cerr << "Read1 after filtering:" << endl;
    finalPostStats->print();

    cerr << endl;
    cerr << "Filtering result:" << endl;
    finalFilterResult->print();

    int *dupHist = NULL;
    double *dupMeanTlen = NULL;
    double *dupMeanGC = NULL;
    double dupRate = 0.0;
    if (mOptions->duplicate.enabled) {
        printf("duplicate enabled is true\n");
        dupHist = new int[mOptions->duplicate.histSize];
        memset(dupHist, 0, sizeof(int) * mOptions->duplicate.histSize);
        dupMeanGC = new double[mOptions->duplicate.histSize];
        memset(dupMeanGC, 0, sizeof(double) * mOptions->duplicate.histSize);
        dupRate = mDuplicate->statAll(dupHist, dupMeanGC, mOptions->duplicate.histSize);
        cerr << endl;
        cerr << "Duplication rate (may be overestimated since this is SE data): " << dupRate * 100.0 << "%" << endl;
    }

    // make JSON report
    JsonReporter jr(mOptions);
    jr.setDupHist(dupHist, dupMeanGC, dupRate);
    jr.report(finalFilterResult, finalPreStats, finalPostStats);

    // make HTML report
    HtmlReporter hr(mOptions);
    hr.setDupHist(dupHist, dupMeanGC, dupRate);
    hr.report(finalFilterResult, finalPreStats, finalPostStats);

    // clean up
    for (int t = 0; t < mOptions->thread; t++) {
        delete threads[t];
        threads[t] = NULL;
        delete configs[t];
        configs[t] = NULL;
    }

    delete finalPreStats;
    delete finalPostStats;
    delete finalFilterResult;

    if (mOptions->duplicate.enabled) {
        delete[] dupHist;
        delete[] dupMeanGC;
    }

    delete[] threads;
    delete[] configs;

    if (leftWriterThread)
        delete leftWriterThread;

    if (!mOptions->split.enabled)
        closeOutput();

    return true;
}

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

    config->totCnt += 1;
    double t0, t1;


    t0 = get_wall_time();

    for (int p = 0; p < pack->count; p++) {

        // original read1
        Read *or1 = pack->data[p];
        // stats the original read before trimming
        t1 = get_wall_time();
        config->getPreStats1()->statRead(or1);//cost 12/57
        config->cost1 += get_wall_time() - t1;


        // handling the duplication profiling
        t1 = get_wall_time();

        //TODO maybe mDuplicate is not thread safe
        if (mDuplicate)
            mDuplicate->statRead(or1);//cost 14/57
        config->cost2 += get_wall_time() - t1;


        t1 = get_wall_time();
        // filter by index
        if (mOptions->indexFilter.enabled && mFilter->filterByIndex(or1)) {
            printf("mOptions->indexFilter ...");
            delete or1;
            continue;
        }
        config->cost3 += get_wall_time() - t1;


        t1 = get_wall_time();
        // umi processing

        //TODO maybe this can be the big hotspot
        if (mOptions->umi.enabled) {
//            printf("mOptions->umi ...");
            mUmiProcessor->process(or1);
        }
        config->cost4 += get_wall_time() - t1;


        t1 = get_wall_time();
        // trim in head and tail, and apply quality cut in sliding window
        //not in
        //TODO look what will do in this function and how to open this function
        //TODO maybe this can be the big hotspot
        Read *r1 = mFilter->trimAndCut(or1, mOptions->trim.front1, mOptions->trim.tail1);
        config->cost5 += get_wall_time() - t1;


        t1 = get_wall_time();
        if (r1 != NULL) {
            //not in because xx is false
            if (mOptions->polyGTrim.enabled) {
                printf("polyGTrim ...\n");
                PolyX::trimPolyG(r1, config->getFilterResult(), mOptions->polyGTrim.minLen);

            }
            //not in because xx is false
            if (mOptions->polyXTrim.enabled) {
                printf("polyXTrim ...\n");
                PolyX::trimPolyX(r1, config->getFilterResult(), mOptions->polyXTrim.minLen);

            }
        }
        config->cost6 += get_wall_time() - t1;


        t1 = get_wall_time();
        //not in because mOptions->adapter.hasSeqR1 is false
        if (r1 != NULL && mOptions->adapter.enabled && mOptions->adapter.hasSeqR1) {
            printf("AdapterTrimmer::trimBySequence ...\n");
            AdapterTrimmer::trimBySequence(r1, config->getFilterResult(), mOptions->adapter.sequence);
        }
        config->cost7 += get_wall_time() - t1;


        t1 = get_wall_time();
        if (r1 != NULL) {
            // not in because mOptions->trim.maxLen1 is 0
            if (mOptions->trim.maxLen1 > 0 && mOptions->trim.maxLen1 < r1->length()) {
                printf("r1->resize ...\n");
                r1->resize(mOptions->trim.maxLen1);
            }
        }
        config->cost8 += get_wall_time() - t1;


        t1 = get_wall_time();
        //in ! O(len)
        int result = mFilter->passFilter(r1);
        config->cost9 += get_wall_time() - t1;


        t1 = get_wall_time();
        config->addFilterResult(result);
        config->cost10 += get_wall_time() - t1;

        //in !
//        if (result != PASS_FILTER) {
//            cout << r1->toString() << endl;
//        }
        if (r1 != NULL && result == PASS_FILTER) {
            t1 = get_wall_time();
            outstr += r1->toString();
            config->cost11 += get_wall_time() - t1;
            // stats the read after filtering
            t1 = get_wall_time();
            config->getPostStats1()->statRead(r1);//cost 23/57
            config->cost12 += get_wall_time() - t1;
            readPassed++;
        }


        t1 = get_wall_time();
        delete or1;
        // if no trimming applied, r1 should be identical to or1
        if (r1 != or1 && r1 != NULL)
            delete r1;
        config->cost13 += get_wall_time() - t1;
    }
    config->cost += get_wall_time() - t0;

    // if splitting output, then no lock is need since different threads write different files
    if (!mOptions->split.enabled)
        mOutputMtx.lock();
    if (mOptions->outputToSTDOUT) {
        printf("mOptions->outputToSTDOUT\n");
        fwrite(outstr.c_str(), 1, outstr.length(), stdout);
    } else if (mOptions->split.enabled) {
        // split output by each worker thread
        if (!mOptions->out1.empty()) {
            printf("mOptions->split.enabled and !mOptions->out1.empty()\n");
            config->getWriter1()->writeString(outstr);
        }
    } else {
        if (mLeftWriter) {
//            printf("mLeftWriter\n");
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

void SingleEndProcessor::initPackRepository() {
    //mRepo.packBuffer = new ReadPack*[PACK_NUM_LIMIT];
//    printf("now initPackRepositorying ...\n");
//    printf("FastqDataChunk totle size %d\n", PACK_NUM_LIMIT);
//    printf("FastqDataChunk pre size %d\n", sizeof(dsrc::fq::FastqDataChunk *));
    mRepo.packBuffer = new dsrc::fq::FastqDataChunk *[PACK_NUM_LIMIT];
    //memset(mRepo.packBuffer, 0, sizeof(ReadPack*)*PACK_NUM_LIMIT);
    memset(mRepo.packBuffer, 0, sizeof(dsrc::fq::FastqDataChunk *) * PACK_NUM_LIMIT);
    mRepo.writePos = 0;
    mRepo.readPos = 0;
    //mRepo.readCounter = 0;

}

void SingleEndProcessor::destroyPackRepository() {
    delete mRepo.packBuffer;
    mRepo.packBuffer = NULL;
}

void SingleEndProcessor::producePack(dsrc::fq::FastqDataChunk *pack) {
    //std::unique_lock<std::mutex> lock(mRepo.mtx);
    /*while(((mRepo.writePos + 1) % PACK_NUM_LIMIT)
        == mRepo.readPos) {
        //mRepo.repoNotFull.wait(lock);
    }*/

    mRepo.packBuffer[mRepo.writePos] = pack;
    //cerr << "chunk in producePack" << endl;
    //cerr << (char *)pack->data.Pointer() << endl;
    mRepo.writePos++;

    /*if (mRepo.writePos == PACK_NUM_LIMIT)
        mRepo.writePos = 0;*/

    //mRepo.repoNotEmpty.notify_all();
    //lock.unlock();
}

void SingleEndProcessor::consumePack(ThreadConfig *config) {
    dsrc::fq::FastqDataChunk *chunk;
    ReadPack *data = new ReadPack;
    //std::unique_lock<std::mutex> lock(mRepo.mtx);
    // buffer is empty, just wait here.
    /*while(mRepo.writePos % PACK_NUM_LIMIT == mRepo.readPos % PACK_NUM_LIMIT) {
        if(mProduceFinished){
            //lock.unlock();
            return;
        }
        //mRepo.repoNotEmpty.wait(lock);
    }*/

    mInputMtx.lock();
    while (mRepo.writePos <= mRepo.readPos) {
        usleep(1000);
        if (mProduceFinished) {
            mInputMtx.unlock();
            return;
        }
    }
    //data = mRepo.packBuffer[mRepo.readPos];
    chunk = mRepo.packBuffer[mRepo.readPos];
    //cerr << "read pos is " << mRepo.readPos << endl;
    //cerr << (char*)chunk->data.Pointer() << endl;

    mRepo.readPos++;

    /*if (mRepo.readPos >= PACK_NUM_LIMIT)
        mRepo.readPos = 0;*/
    mInputMtx.unlock();


    double t = get_wall_time();
    //data format for from dsrc to fastp
    data->count = dsrc::fq::chunkFormat(chunk, data->data, true);

    config->costFormat += get_wall_time() - t;
    //cerr << (char*)chunk->data.Pointer() << endl;
    fastqPool->Release(chunk);

    //lock.unlock();
    //mRepo.repoNotFull.notify_all();

    processSingleEnd(data, config);


}

void SingleEndProcessor::producerTask() {
    if (mOptions->verbose)
        loginfo("start to load data");
    long lastReported = 0;
    int slept = 0;
    long readNum = 0;
    bool splitSizeReEvaluated = false;
    //Read** data = new Read*[PACK_SIZE];
    //memset(data, 0, sizeof(Read*)*PACK_SIZE);
    //FastqReader reader(mOptions->in1, true, mOptions->phred64);
    //dsrc::fq::FastqDataPool* fastqPool = new dsrc::fq::FastqDataPool(32,1<<22);
    dsrc::fq::FastqFileReader *fileReader = new dsrc::fq::FastqFileReader(mOptions->in1);
    dsrc::fq::FastqReader *dataReader = new dsrc::fq::FastqReader(*fileReader, *fastqPool);


    dsrc::fq::FastqDataChunk *chunk;
    while ((chunk = dataReader->readNextChunk()) != NULL) {
        //cerr << "chunk content in producer =======================" << endl;
        //cerr << (char*) chunk->data.Pointer() << endl;
        producePack(chunk);
        //cerr << (char*)mRepo.packBuffer[0]->data.Pointer() << endl;
        //usleep(200);
        while (mRepo.writePos - mRepo.readPos > PACK_IN_MEM_LIMIT) {
            //cerr<<"sleep"<<endl;
            slept++;
            usleep(100);
        }

    }



    //std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
    mProduceFinished = true;
    if (mOptions->verbose)
        loginfo("all reads loaded, start to monitor thread status");
    //lock.unlock();
    delete fileReader;
    delete dataReader;
}

void SingleEndProcessor::consumerTask(ThreadConfig *config) {
//    cost1 = 0;
//    cost2 = 0;
//    cost3 = 0;
//    cost4 = 0;
//    cost5 = 0;
//    cost6 = 0;
//    cost7 = 0;
//    cost8 = 0;
//    cost9 = 0;
//    cost10 = 0;
//    cost11 = 0;
//    cost11b = 0;
//    cost12 = 0;
//    totCnt = 0;
    //writePos is the postion where producer has written in buffer
    //readPos is the posthon where consumers has read from buffer
    //if mRepo.writePos <= mRepo.readPos, whether all the tasks has down or consumers have to wait producer.

    while (true) {
        if (config->canBeStopped()) {
            mFinishedThreads++;
            break;
        }
        while (mRepo.writePos <= mRepo.readPos) {
            if (mProduceFinished)
                break;
            usleep(1000);
        }
        //std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
        if (mProduceFinished && mRepo.writePos == mRepo.readPos) {
            mFinishedThreads++;
            if (mOptions->verbose) {
                string msg = "thread " + to_string(config->getThreadId() + 1) + " data processing completed";
                loginfo(msg);
            }
            //lock.unlock();
            break;
        }

//        assert(mRepo.writePos > mRepo.readPos);
        if (mOptions->verbose) {
            string msg = "thread " + to_string(config->getThreadId() + 1) + " is processing the " +
                         to_string(mRepo.readPos) + " / " + to_string(mRepo.writePos) + " pack";
            loginfo(msg);
        }
        consumePack(config);

        //TODO if ?? maybe 3GS
//        if (mProduceFinished) {
//            if (mOptions->verbose) {
//                string msg = "thread " + to_string(config->getThreadId() + 1) + " is processing the " +
//                             to_string(mRepo.readPos) + " / " + to_string(mRepo.writePos) + " pack";
//                loginfo(msg);
//            }
//            consumePack(config);
//            //lock.unlock();
//        } else {
//            //lock.unlock();
//            consumePack(config);
//        } //--[haoz:] I think it 3GS can add here
    }
//
//
//    printf("total cost1 : %.5f\n", cost1);
//    printf("total cost2 : %.5f\n", cost2);
//    printf("total cost3 : %.5f\n", cost3);
//    printf("total cost4 : %.5f\n", cost4);
//    printf("total cost5 : %.5f\n", cost5);
//    printf("total cost6 : %.5f\n", cost4);
//    printf("total cost7 : %.5f\n", cost7);
//    printf("total cost8 : %.5f\n", cost8);
//    printf("total cost9 : %.5f\n", cost9);
//    printf("total cost10 : %.5f\n", cost10);
//    printf("total cost11b : %.5f\n", cost11b);
//    printf("total cost11 : %.5f\n", cost11);
//    printf("total cost12 : %.5f\n", cost12);
//    printf("total statReadCost : %.5f\n", cost11 + cost11);
//    printf("total costTotel : %.5f\n",
//           cost1 + cost2 + cost3 + cost4 + cost5 + cost6 + cost7 + cost8 + cost9 + cost10 + cost11b + cost11 + cost12);
//    printf("total cost : %.5f\n", cost);
//    printf("total  : %d\n", totCnt);
//
//
//
////TODO check another !
//    fstream out;
//    fstream in;
//
//    //TODO it seems that the output of file1 is different everytime,which only happen when read data with adapter.
//    string outFileName1 = "oneStatsCheckOfThread" + to_string(config->getThreadId());
//    string outFileName2 = "kMerCheckOfThread" + to_string(config->getThreadId());
//    long *tmpKmer = config->getPreStats1()->mKmer;
//    int tmpKmerSize = config->getPreStats1()->mKmerBufLen;
//    long *tmpStats = config->getPreStats1()->mCycleTotalQual;
//    int tmpStatsSize = config->getPreStats1()->mBufLen;
//    mCounts = new uint16[mKeyLenInBit];
//    uint16 tmpMcounts=
//    printf("Start checking ...");
//    printf("Open file : %s\n", outFileName1.c_str());
//
//    // Read sampling file
//    long *temp = new long[tmpStatsSize];
//    in.open(outFileName1.c_str(), ios::in | ios::binary);
//    if (!in) {
//        printf("Can't open file \"%s\"\n", outFileName1.c_str());
//    } else {
//        in.seekg(0, ios::beg);
//        in.read(reinterpret_cast<char *>(temp), tmpStatsSize * sizeof(long));
//        printf("=================================================\n");
//        for (int i = 0; i < tmpStatsSize; i++) {
//            if (temp[i] != tmpStats[i]) {
//                printf("GG on test %d  STD : %ld   Now : %ld\n", i, temp[i], tmpStats[i]);
//            }
////            printf("%ld\n", temp[i]);
//        }
//        printf("=================================================\n");
//    }
//
//
//    printf("=================================================\n");
//    for (int i = 0; i < tmpStatsSize; i++) {
//        printf("%ld ", tmpStats[i]);
//        if ((i + 1) % 20 == 0)printf("\n");
//    }
//    printf("=================================================\n");
//
//    out.open(outFileName1.c_str(), ios::out | ios::binary);
//    out.seekp(0, ios::beg);
//    out.write(reinterpret_cast<char *>(tmpStats), tmpStatsSize * sizeof(long));
//    out.close();
//
//    out.open(outFileName2.c_str(), ios::out | ios::binary);
//    out.seekp(0, ios::beg);
//    out.write(reinterpret_cast<char *>(tmpKmer), tmpKmerSize * sizeof(long));
//    out.close();


    if (mFinishedThreads == mOptions->thread) {
        if (mLeftWriter)
            mLeftWriter->setInputCompleted();
    }

    if (mOptions->verbose) {
        string msg = "thread " + to_string(config->getThreadId() + 1) + " finished";
        loginfo(msg);
    }

}

void SingleEndProcessor::writeTask(WriterThread *config) {
    while (true) {
        if (config->isCompleted()) {
            // last check for possible threading related issue
            config->output();
            break;
        }
        config->output();
    }

    if (mOptions->verbose) {
        string msg = config->getFilename() + " writer finished";
        loginfo(msg);
    }
}
