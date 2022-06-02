#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <random>
#include <unordered_map>
#include <algorithm>
#include <ctime>
#include "MurmurHash3.h"

using namespace std;

#define KEY_SIZE 13

//the exp part and coe part share 16 bits, and bot of them cannot be less than 1.
#define EXP_PART_BITS 2
#define EXP_PART_MASK 3

//these are corresponding definitions to reduce calculation
#define COE_PART_BITS 14
#define COE_PART_MASK 0x3FFF
#define COE_PART_BASEVAL 0x2000

#define DECAY_PARAM 1.001

#define HASH_SEED 0

static mt19937 rng(time(0));

bool cmpPairFunc(pair<uint8_t *, int> p1, pair<uint8_t *, int> p2) {
    return p1.second > p2.second;
}

//for recording information in the unordered_map
struct CmpFunc {
    bool operator()(const uint8_t *keyA, const uint8_t *keyB) const {
        return memcmp(keyA, keyB, KEY_SIZE) == 0;
    }
};

struct HashFunc {
    unsigned int operator()(const uint8_t* key) const {
        unsigned int hashValue=0;
        MurmurHash3_x86_32(key, KEY_SIZE,0, &hashValue);;
        return hashValue;
    }
};

struct HeapNode {
    uint8_t key[KEY_SIZE]={0};
    uint32_t counter=0;
};

template<uint32_t HEAP_MAX_SIZE>
class MinHeap {
private:
    HeapNode array[HEAP_MAX_SIZE];
    uint32_t maxVal = 0;//record the max value in minHeap
    uint32_t size = 0;//record the current size of minHeap
public:
    MinHeap() {}

    uint32_t getMinVal() {
        if (size > 0) {
            return array[0].counter;
        } else {
            return -1;
        }
    }

    uint32_t getMaxVal() {
        return maxVal;
    }

    uint32_t getSize(){
        return size;
    }

    uint8_t* getKey(uint32_t nodeIndex){
        return array[nodeIndex].key;
    }
    uint32_t getValue(uint32_t nodeIndex){
        return array[nodeIndex].counter;
    }

    bool findKey(uint8_t *key, uint32_t &nodeIndex, uint32_t &counter) {
        for (uint32_t i = 0; i < size; i++) {
            if (memcmp(array[i].key, key, KEY_SIZE) == 0) {
                nodeIndex = i;
                counter = array[i].counter;
                return true;
            }
        }
        return false;
    }

    void heapPush(uint8_t *key, uint32_t counter) {
        if (size < HEAP_MAX_SIZE) {
            memcpy(array[size].key, key, KEY_SIZE);
            array[size].counter = counter;
            heapRepairUp(size);
            size++;

            maxVal = max(maxVal, counter);
        } else {//heap is full
            if (array[0].counter < counter) {
                memcpy(array[0].key, key, KEY_SIZE);
                array[0].counter = counter;
                heapRepairDown(0);

                maxVal = max(maxVal, counter);
            }
        }
    }

    void updateHeapNode(uint32_t nodeIndex,uint32_t newCounter) {//only consider the case that new counter is larger than old counter
        array[nodeIndex].counter = newCounter;
        heapRepairDown(nodeIndex);
        maxVal = max(maxVal, newCounter);
    }

private:
    void heapRepairDown(uint32_t nodeIndex) {
        uint32_t l, r;
        uint32_t smallestNodeIdx;

        while (nodeIndex < (size >> 1)) {
            l = (nodeIndex << 1) + 1;
            r = (nodeIndex + 1) << 1;

            smallestNodeIdx = nodeIndex;
            if (l < size && array[l].counter < array[smallestNodeIdx].counter) {
                smallestNodeIdx = l;
            }
            if (r < size && array[r].counter < array[smallestNodeIdx].counter) {
                smallestNodeIdx = r;
            }

            if (smallestNodeIdx != nodeIndex) {
                HeapNode temp = array[nodeIndex];
                array[nodeIndex] = array[smallestNodeIdx];
                array[smallestNodeIdx] = temp;

                nodeIndex = smallestNodeIdx;
            } else {
                break;
            }
        }
    }

    void heapRepairUp(uint32_t nodeIndex) {
        uint32_t parentIdx;
        while (nodeIndex > 0) {
            parentIdx = (nodeIndex - 1) >> 1;

            if (array[nodeIndex].counter < array[parentIdx].counter) {
                HeapNode temp = array[nodeIndex];
                array[nodeIndex] = array[parentIdx];
                array[parentIdx] = temp;
                nodeIndex = parentIdx;
            } else {
                break;
            }
        }
    }
};

//get the real value represented by the two mode active counter
inline uint32_t getRealValOfTwoMoldActiveCounter(uint16_t counter) {
    if (counter > 0x8000) {//in active mode
        uint32_t coePart = counter >> EXP_PART_BITS;
        uint32_t expPart = counter & EXP_PART_MASK;
        uint32_t value = coePart << (EXP_PART_BITS + expPart);
        return value;
    } else {//in normal mode
        return counter;
    }
}

inline uint16_t incTwoMoldActiveCounter(uint16_t counter) {
    if (counter >= 0x8000) {//in active mode
        uint32_t coePart = counter >> EXP_PART_BITS;
        uint32_t expPart = counter & EXP_PART_MASK;

        double randomVal = (double)rng()/4294967296;
        if (randomVal * ((uint32_t) 1 << (EXP_PART_BITS + expPart)) < 1.0) {
            if (coePart < COE_PART_MASK) {//the coe part can be directly increased by one without overflow
                coePart += 1;
            } else {
                if (expPart < EXP_PART_MASK) {
                    coePart = COE_PART_BASEVAL;
                    expPart += 1;
                }
            }
            uint16_t newCounter = (coePart << EXP_PART_BITS) + expPart;
            return newCounter;
        } else {
            return counter;
        }
    } else {//in normal mode
        return counter + 1;
    }
}

inline uint16_t addTwoMoldActiveCounter(uint16_t counter, uint32_t &realVal) {
    if (counter >= 0x8000) {//in active mode
        uint32_t coePart = counter >> EXP_PART_BITS;
        uint32_t expPart = counter & EXP_PART_MASK;

        if (coePart < COE_PART_MASK) {//the coe part can be directly increased by one without overflow
            coePart += 1;
        } else {
            if (expPart < EXP_PART_MASK) {
                coePart = COE_PART_BASEVAL;
                expPart += 1;
            }
        }
        realVal = coePart << (EXP_PART_BITS + expPart);
        uint16_t newCounter = (coePart << EXP_PART_BITS) + expPart;
        return newCounter;
    } else {//in normal mode
        realVal = counter + 1;
        return realVal;
    }
}

inline uint16_t decTwoMoldActiveCounter(uint16_t counter) {
    if (counter > 0x8000) {//after decreasing, the counter is still in active mode
        uint32_t coePart = counter >> EXP_PART_BITS;
        uint32_t expPart = counter & EXP_PART_MASK;

        if (coePart > COE_PART_BASEVAL) {
            coePart -= 1;
        } else {
            coePart = COE_PART_MASK;
            expPart -= 1;
        }

        double randomVal = (double)rng()/4294967296;
        if (randomVal * ((uint32_t) 1 << (EXP_PART_BITS + expPart)) < 1.0) {
            uint16_t newCounter = (coePart << EXP_PART_BITS) + expPart;
            return newCounter;
        } else {
            return counter;
        }
    } else {//in normal mode
        return counter - 1;
    }
}

struct Bucket {
    uint16_t fingerprint = 0;
    uint16_t activeCounter = 0;
};

//ROW_NUM cannot be larger than 3
template<uint32_t ROW_NUM, uint32_t COL_NUM, uint32_t HEAP_MAX_SIZE>
class ActiveKeeper {
private:
    Bucket bucketArray[ROW_NUM][COL_NUM];
    MinHeap<HEAP_MAX_SIZE> minHeap;
public:
    ActiveKeeper() : minHeap() {}

    void insert(uint8_t *key) {
        uint32_t findedNodeIndex;
        uint32_t counterInMinHeap;
        bool flag;
        flag = minHeap.findKey(key, findedNodeIndex, counterInMinHeap);

        if (flag) {
            minHeap.updateHeapNode(findedNodeIndex, counterInMinHeap + 1);
            counterInMinHeap++;
        }

        uint64_t hashValue[2]{0};
        MurmurHash3_x64_128(key, KEY_SIZE, HASH_SEED, &hashValue);

        uint16_t fingerprint = hashValue[1] & 0xffff;
        uint32_t estimatedVal = 0;

        uint32_t *hashValues = (uint32_t *) (&hashValue);
        uint32_t minValInMinHeap = minHeap.getMinVal();
        uint32_t maxValInMinHeap = minHeap.getMaxVal();

        for (uint32_t rowIndex = 0; rowIndex < ROW_NUM; rowIndex++) {
            uint32_t colIndex = hashValues[rowIndex] % COL_NUM;

            Bucket &bucket = bucketArray[rowIndex][colIndex];
            if (bucket.activeCounter == 0) {
                bucket.fingerprint = fingerprint;
                bucket.activeCounter = 1;
                estimatedVal = max(estimatedVal, (uint32_t) 1);
            } else if (bucket.fingerprint == fingerprint) {
                uint32_t ACCounterVal = getRealValOfTwoMoldActiveCounter(bucket.activeCounter);
                uint32_t addedACCounterVal;
                uint16_t addedACCounter = addTwoMoldActiveCounter(bucket.activeCounter, addedACCounterVal);

                if (flag && addedACCounterVal <= counterInMinHeap) {
                    bucket.activeCounter = addedACCounter;
                } else if (!flag && ACCounterVal <= minValInMinHeap) {
                    if (addedACCounterVal > minValInMinHeap + 1) {
                        double randomVal = (double)rng()/4294967296;
                        if (randomVal < (1.0 / (minValInMinHeap + 1 - ACCounterVal))) {
                            estimatedVal = max(estimatedVal, minValInMinHeap + 1);
                        } else {
                            estimatedVal = max(estimatedVal, ACCounterVal);
                        }
                    } else {
                        bucket.activeCounter = incTwoMoldActiveCounter(bucket.activeCounter);
                        ACCounterVal = getRealValOfTwoMoldActiveCounter(bucket.activeCounter);
                        estimatedVal = max(estimatedVal, ACCounterVal);
                    }
                }
            } else {
                uint32_t ACCounterVal = getRealValOfTwoMoldActiveCounter(bucket.activeCounter);

                double randomVal = (double)rng()/4294967296;
                double prob=pow(DECAY_PARAM +(double) minValInMinHeap /maxValInMinHeap,-1.0 *ACCounterVal);
                if ((ACCounterVal <= minValInMinHeap || !flag) && randomVal <prob){
                    uint16_t newACCounter = decTwoMoldActiveCounter(bucket.activeCounter);
                    if (newACCounter == 0) {
                        bucket.fingerprint = fingerprint;
                        bucket.activeCounter = 1;
                        estimatedVal = max(estimatedVal, (uint32_t) 1);
                    } else {
                        bucket.activeCounter = newACCounter;
                    }
                }
            }
        }

        if (!flag) {
            minHeap.heapPush(key, estimatedVal);
        }
    }

    void getEstimatedFlowSizes(unordered_map<uint8_t *, unsigned int, HashFunc, CmpFunc> &estimatedFlowSizes) {
        for (uint32_t i = 0; i < minHeap.getSize(); i++) {
            uint8_t *tempKey = (uint8_t *) malloc(KEY_SIZE);
            memcpy(tempKey, minHeap.getKey(i), KEY_SIZE);
            estimatedFlowSizes[tempKey] = minHeap.getValue(i);
        }
    }
};

unsigned int ReadInTraces(const char* tracePreFix, uint8_t** keys, unordered_map<uint8_t*, unsigned int, HashFunc, CmpFunc>& actualFlowSizes, unsigned int maxItemNum)
{
    unsigned int count = 0;

    unsigned int countInFile = 0;
    for (int datafileCnt = 0; datafileCnt <= 10; ++datafileCnt){
        countInFile = 0;
        char traceFilePath[100];
        sprintf(traceFilePath, "%s%d.dat", tracePreFix, datafileCnt);
        printf("Start reading %s\n", traceFilePath);

        FILE* fin = fopen(traceFilePath, "rb");
        char temp[KEY_SIZE]{ 0 };
        uint8_t* key;
        while (fread(temp, 1, KEY_SIZE, fin) == KEY_SIZE) {
            key = (uint8_t*)malloc(KEY_SIZE);
            memcpy(key, temp, KEY_SIZE);
            keys[count] = key;
            if (actualFlowSizes.find(key) == actualFlowSizes.end()) {
                actualFlowSizes[key] = 1;
            }
            else {
                actualFlowSizes[key] += 1;
            }
            count++;
            if(count>=maxItemNum){
                printf("The dataset has more than %d items, set a larger value for maxItemNum", maxItemNum);
                exit(-1);
            }

            countInFile++;
            if (countInFile % 5000000 == 0) {
                printf("\thave read %u items in %s, the dataset now has %u items\n", countInFile,traceFilePath, count);
            }
        }
        fclose(fin);
        printf("Finish reading %s (%u items), the dataset now has %u items\n", traceFilePath,countInFile,count);

    }
	return count;
}



int main()
{
    //prepare dataset
    cout << "prepare dataset" << endl;
    unsigned int maxItemNum = 40*1000000;//max number of items
    uint8_t** keys = (uint8_t**)calloc(maxItemNum, sizeof(uint8_t*));//to store keys of all items
    unordered_map<uint8_t*, unsigned int, HashFunc, CmpFunc> actualFlowSizes;//count the actual flow sizes
    unsigned int actualItemNum = ReadInTraces(R"(../../data/)", keys, actualFlowSizes,maxItemNum);
    cout << "number of items: " << actualItemNum << endl;
    cout << "number of flows: " << actualFlowSizes.size() << endl;
    cout << "*********************" << endl;

    //prepare algorithm
    cout << "prepare algorithm" << endl;

    //parameters

    const uint32_t rowNum = 2;
    //top-1000
    const uint32_t heapMaxSize = 1000;
    const uint32_t colNum = 23694;//200KB

    ActiveKeeper<rowNum, colNum, heapMaxSize> activeKeeper;

    //output memory info
    double bucketMem = 4 * rowNum * colNum / 1024.0;//32 bits for each bucket
    double minHeapMem = heapMaxSize * (KEY_SIZE + 18 / 8.0) / 1024.0;//122 bits for each slot
    cout << "bucketMem: " << bucketMem << "KB" << endl;
    cout << "minHeapMem: " << minHeapMem << "KB" << endl;
    cout << "totalMem:" << bucketMem + minHeapMem << "KB" << endl;
    cout << "*********************" << endl;

    //insert items
    cout << "insert items" << endl;
    clock_t time1 = clock();
    for (unsigned int i = 0; i < actualItemNum; i++) {
        activeKeeper.insert(keys[i]);
    }
    clock_t time2 = clock();

    //calculate throughput
    double numOfSeconds = (double) (time2 - time1) / CLOCKS_PER_SEC;//the seconds using to insert items
    double throughput = (actualItemNum / 1000000.0) / numOfSeconds;
    cout << "use " << numOfSeconds << " seconds" << endl;
    cout << "throughput: " << throughput << " Mpps" << ", each insert operation uses " << 1000.0 / throughput << " ns"
         << endl;
    cout << "*********************" << endl;

    //calculate precision, ARE, AAE
    //get sorted acutal flow sizes and sorted estimated flow sizes
    vector<pair<uint8_t *, unsigned int>> actualFlowSizesVector;
    for (auto iter = actualFlowSizes.begin(); iter != actualFlowSizes.end(); iter++) {
        actualFlowSizesVector.push_back(make_pair(iter->first, iter->second));
    }
    sort(actualFlowSizesVector.begin(), actualFlowSizesVector.end(), cmpPairFunc);

    vector<pair<uint8_t *, unsigned int>> estimatedFlowSizesVector;
    unordered_map<uint8_t *, unsigned int, HashFunc, CmpFunc> estimatedFlowSizes;
    activeKeeper.getEstimatedFlowSizes(estimatedFlowSizes);

    for (auto iter = estimatedFlowSizes.begin(); iter != estimatedFlowSizes.end(); iter++) {
        estimatedFlowSizesVector.push_back(make_pair(iter->first, iter->second));
    }
    cout << "get top-" << estimatedFlowSizes.size() << " flows" << endl;
    sort(estimatedFlowSizesVector.begin(), estimatedFlowSizesVector.end(), cmpPairFunc);

    uint32_t k = min((unsigned int) heapMaxSize, (unsigned int) estimatedFlowSizes.size());
    //get acutal top-k flows
    unordered_map<uint8_t *, unsigned int, HashFunc, CmpFunc> actualTopKFlowSizes;
    for (unsigned int i = 0; i < k; i++) {
        unsigned int actualSize = actualFlowSizesVector[i].second;
        actualTopKFlowSizes[actualFlowSizesVector[i].first] = actualSize;
    }

    cout << "Top-" << k << " flow detection Metrics" << endl;
    k = min((unsigned int) estimatedFlowSizes.size(), k);
    //top-k precision
    unsigned int TPNumOfTopK = 0;//True positive num
    for (int i = 0; i < k; i++) {
        uint8_t *key = estimatedFlowSizesVector[i].first;
        if (actualTopKFlowSizes.find(key) != actualTopKFlowSizes.end()) {
            TPNumOfTopK++;
        }
    }
    cout << "top-"<<k<<" Precision: " << (double) TPNumOfTopK / k << endl;

    //ARE,AAE of the reported top-k flows
    double totalREOfTopK = 0;//total relative error of top-k flows
    double totalAEOfTopK = 0;
    for (int i = 0; i < k; i++) {
        uint8_t *key = estimatedFlowSizesVector[i].first;
        unsigned int estimatedSize = estimatedFlowSizesVector[i].second;
        unsigned int actualSize = 0;
        if (actualFlowSizes.find(key) != actualFlowSizes.end()) {
            actualSize = actualFlowSizes[key];
        }

        totalAEOfTopK += abs((double) estimatedSize - actualSize);
        double rE = abs((double) estimatedSize - actualSize) / actualSize;
        totalREOfTopK += rE;
    }
    cout << "ARE of reported top-k flows: " << totalREOfTopK / k << endl;
    cout << "AAE of reported top-k flows: " << totalAEOfTopK / k << endl;

    //release resources
    for (auto iter = estimatedFlowSizes.begin(); iter != estimatedFlowSizes.end(); iter++) {
        free(iter->first);
    }
    for (unsigned int i = 0; i < actualItemNum; i++) {
        free(keys[i]);
    }
    free(keys);
    return 0;
}