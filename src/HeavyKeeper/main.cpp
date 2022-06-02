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
#define DECAY_PARAM 1.08

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

struct SSKeyNode;

struct SSValueNode {
    uint32_t counter;
    SSValueNode *next;
    SSValueNode *pre;
    SSKeyNode *firstKeyNode;
};

struct SSKeyNode {
    uint8_t key[KEY_SIZE];
    SSValueNode *parent;
    SSKeyNode *next;
};

template<uint32_t SS_MAX_SIZE>
class StreamSummary {
private:
    uint32_t size = 0;
    SSValueNode *firstValueNode = nullptr;
    unordered_map<uint8_t *, SSKeyNode *, HashFunc, CmpFunc> hashTable;
public:
    StreamSummary() {}
    ~StreamSummary(){
        SSValueNode *ssValueNode = firstValueNode;
        while (ssValueNode != nullptr) {
            SSKeyNode *ssKeyNode = ssValueNode->firstKeyNode;
            do {
                SSKeyNode *tempKeyNode = ssKeyNode;
                ssKeyNode = ssKeyNode->next;
                free(tempKeyNode);
            } while (ssKeyNode != ssValueNode->firstKeyNode);

            SSValueNode *tempValueNode = ssValueNode;
            ssValueNode = ssValueNode->next;
            free(tempValueNode);
        }
    }

    uint32_t getSize() {
        return size;
    }
    uint32_t getMinVal() {
        if (size != 0) {
            return firstValueNode->counter;
        } else {
            return -1;
        }
    }

    SSValueNode* getFirstValueNode(){
        return firstValueNode;
    }

    SSKeyNode *findKey(uint8_t *key) {
        auto iter = hashTable.find(key);
        if (iter == hashTable.end()) {
            return nullptr;
        } else {
            return iter->second;
        }
    }


    void SSPush(uint8_t *key, uint32_t counter) {
        if (size == 0) {
            SSValueNode *ssValueNode = (SSValueNode *) calloc(1, sizeof(SSValueNode));

            SSKeyNode *ssKeyNode = (SSKeyNode *) calloc(1, sizeof(SSKeyNode));
            memcpy(ssKeyNode->key, key, KEY_SIZE);
            ssKeyNode->parent = ssValueNode;
            ssKeyNode->next = ssKeyNode;

            ssValueNode->counter = counter;
            ssValueNode->firstKeyNode = ssKeyNode;

            firstValueNode = ssValueNode;

            hashTable[ssKeyNode->key] = ssKeyNode;
            size++;
        } else if (size < SS_MAX_SIZE) {//in this case of heavy keeper, the new node will not be larger than minVal
            if (counter == firstValueNode->counter) {//add to the first value node, behind the first key node
                SSKeyNode *ssKeyNode = (SSKeyNode *) calloc(1, sizeof(SSKeyNode));
                memcpy(ssKeyNode->key, key, KEY_SIZE);
                ssKeyNode->parent = firstValueNode;
                ssKeyNode->next = firstValueNode->firstKeyNode->next;

                firstValueNode->firstKeyNode->next = ssKeyNode;

                hashTable[ssKeyNode->key] = ssKeyNode;
                size++;
            } else {//create a new value node and add it before the first value node
                SSValueNode *ssValueNode = (SSValueNode *) calloc(1, sizeof(SSValueNode));

                SSKeyNode *ssKeyNode = (SSKeyNode *) calloc(1, sizeof(SSKeyNode));
                memcpy(ssKeyNode->key, key, KEY_SIZE);
                ssKeyNode->parent = ssValueNode;
                ssKeyNode->next = ssKeyNode;

                ssValueNode->counter = counter;
                ssValueNode->next = firstValueNode;
                ssValueNode->firstKeyNode = ssKeyNode;

                firstValueNode->pre = ssValueNode;
                firstValueNode = ssValueNode;

                hashTable[ssKeyNode->key] = ssKeyNode;
                size++;
            }
        } else if (counter >
                   firstValueNode->counter) {//when stream summary is full, in heavy keeper, the new pushed value must be minVal+1
            SSKeyNode *ssKeyNode = firstValueNode->firstKeyNode;

            if (ssKeyNode->next == ssKeyNode) {//the first value node only has one key node
                //delete the first key node
                auto iter = hashTable.find(ssKeyNode->key);
                hashTable.erase(iter);
                memcpy(ssKeyNode->key, key, KEY_SIZE);

                SSValueNode *nextValueNode = firstValueNode->next;
                if (nextValueNode != nullptr && counter ==
                                                nextValueNode->counter) {//move the key node to the next value node, and delete the old value node
                    free(firstValueNode);

                    ssKeyNode->parent = nextValueNode;
                    ssKeyNode->next = nextValueNode->firstKeyNode->next;
                    nextValueNode->firstKeyNode->next = ssKeyNode;

                    firstValueNode = nextValueNode;
                    hashTable[ssKeyNode->key] = ssKeyNode;
                } else {
                    firstValueNode->counter = counter;
                    hashTable[ssKeyNode->key] = ssKeyNode;
                }
            } else {
                ssKeyNode = firstValueNode->firstKeyNode->next;
                //delete the second key node
                auto iter = hashTable.find(ssKeyNode->key);
                hashTable.erase(iter);
                memcpy(ssKeyNode->key, key, KEY_SIZE);

                firstValueNode->firstKeyNode->next = ssKeyNode->next;

                SSValueNode *nextValueNode = firstValueNode->next;
                if (nextValueNode != nullptr &&
                    counter == nextValueNode->counter) {//move the key node to the next value node
                    ssKeyNode->parent = nextValueNode;
                    ssKeyNode->next = nextValueNode->firstKeyNode->next;
                    nextValueNode->firstKeyNode->next = ssKeyNode;

                    hashTable[ssKeyNode->key] = ssKeyNode;
                } else {//create a new value node
                    SSValueNode *ssValueNode = (SSValueNode *) calloc(1, sizeof(SSValueNode));

                    ssKeyNode->parent = ssValueNode;
                    ssKeyNode->next = ssKeyNode;

                    ssValueNode->counter = counter;
                    ssValueNode->next = nextValueNode;
                    ssValueNode->pre = firstValueNode;
                    ssValueNode->firstKeyNode = ssKeyNode;

                    firstValueNode->next = ssValueNode;
                    if (nextValueNode != nullptr) {
                        nextValueNode->pre = ssValueNode;
                    }

                    hashTable[ssKeyNode->key] = ssKeyNode;
                }
            }
        }
    }

//in HeavyKeeper, the newCounter may be larger than oldCounter+1, though the situation is rare.
//for example: Suppose f1 and f2 share the same bucket and fingerprint. When the counter equal to minVal+1, f1 is added into stream summary.
// Suppose the counter is decayed, and f2 increase this counter to minVal+1 again. This time, f2 is added into stream summary.
// Afterwards, the update of f1 or f2 is the sum of their increments.
    void SSUpdate(SSKeyNode *ssKeyNode, uint32_t newCounter) {
        SSValueNode *ssValueNode = ssKeyNode->parent;

        if (ssKeyNode->next == ssKeyNode) {//the value node only has one key node
            SSValueNode *candidateValueNode = ssValueNode->next;
            while (candidateValueNode != nullptr && candidateValueNode->counter <
                                                    newCounter) {//find the first value node that is not smaller than newCounter
                candidateValueNode = candidateValueNode->next;
            }

            if (candidateValueNode != nullptr && candidateValueNode->counter ==
                                                 newCounter) {//find the value node that equals to newCounter, need to remove the old value node and move the key node
                //remove old valueNode
                if (firstValueNode == ssValueNode) {
                    firstValueNode = ssValueNode->next;
                } else {
                    ssValueNode->next->pre = ssValueNode->pre;
                    ssValueNode->pre->next = ssValueNode->next;
                }
                free(ssValueNode);

                //move the key node
                ssKeyNode->parent = candidateValueNode;
                ssKeyNode->next = candidateValueNode->firstKeyNode->next;
                candidateValueNode->firstKeyNode->next = ssKeyNode;
            } else {//do not find the equal node, just need to modify the current value node
                ssValueNode->counter = newCounter;
            }
        } else {//the value node has other key nodes, thus the old value node cannot be removed
            //detach the key node from old value node
            //Since the linked list is single linked, we have to traverse the list to find the pre node.
            //We choose to replace the node content instead of finding the pre node.
            SSKeyNode *nextKeyNode = ssKeyNode->next;

            auto iter = hashTable.find(ssKeyNode->key);
            hashTable.erase(iter);
            iter = hashTable.find(nextKeyNode->key);
            hashTable.erase(iter);

            uint8_t key[KEY_SIZE];
            memcpy(key, ssKeyNode->key, KEY_SIZE);
            memcpy(ssKeyNode->key, nextKeyNode->key, KEY_SIZE);
            memcpy(nextKeyNode->key, key, KEY_SIZE);

            hashTable[ssKeyNode->key] = ssKeyNode;
            hashTable[nextKeyNode->key] = nextKeyNode;

            ssKeyNode->next = nextKeyNode->next;
            if (ssValueNode->firstKeyNode == nextKeyNode) {
                ssValueNode->firstKeyNode = ssKeyNode;
            }
            ssKeyNode = nextKeyNode;


            SSValueNode *preOfCandidateValueNode = ssValueNode;
            SSValueNode *candidateValueNode = ssValueNode->next;

            while (candidateValueNode != nullptr && candidateValueNode->counter < newCounter) {
                preOfCandidateValueNode = candidateValueNode;
                candidateValueNode = candidateValueNode->next;
            }

            if (candidateValueNode != nullptr && candidateValueNode->counter == newCounter) {
                //move the key node
                ssKeyNode->parent = candidateValueNode;
                ssKeyNode->next = candidateValueNode->firstKeyNode->next;
                candidateValueNode->firstKeyNode->next = ssKeyNode;
            } else {
                SSValueNode *newValueNode = (SSValueNode *) calloc(1, sizeof(SSValueNode));
                ssKeyNode->parent = newValueNode;
                ssKeyNode->next = ssKeyNode;

                newValueNode->counter = newCounter;
                newValueNode->firstKeyNode = ssKeyNode;
                newValueNode->next = candidateValueNode;
                newValueNode->pre = preOfCandidateValueNode;

                preOfCandidateValueNode->next = newValueNode;
                if (candidateValueNode != nullptr) {
                    candidateValueNode->pre = newValueNode;
                }
            }
        }
    }
};


struct Bucket {
    uint16_t fingerprint = 0;
    uint32_t counter = 0;
};

//minimum version
//ROW_NUM cannot be larger than 4
template<uint32_t ROW_NUM, uint32_t COL_NUM, uint32_t SS_MAX_SIZE>
class HeavyKeeper {
private:
    Bucket bucketArray[ROW_NUM][COL_NUM];
    StreamSummary<SS_MAX_SIZE> ss;
public:
    HeavyKeeper() {};

    inline void insert(uint8_t *key) {
        SSKeyNode *findedKeyNodePtr;
        findedKeyNodePtr = ss.findKey(key);

        uint64_t hashValue[2]{0};
        MurmurHash3_x64_128(key, KEY_SIZE, HASH_SEED,&hashValue);

        uint16_t fingerprint = hashValue[1] & 0xffff;

        uint32_t *hashValues = (uint32_t *) (&hashValue);
        uint32_t minValInSS = ss.getMinVal();
        Bucket *firstEmptyBucketPtr = nullptr;
        uint32_t counterOfMinBucket = -1;
        Bucket *minBucketPtr = nullptr;

        bool addFlag = false;
        uint32_t toUpdateValue = 0;

        for (uint32_t rowIndex = 0; rowIndex < ROW_NUM; rowIndex++) {
            uint32_t colIndex = hashValues[rowIndex] % COL_NUM;

            Bucket &bucket = bucketArray[rowIndex][colIndex];

            if (bucket.counter == 0) {
                if (firstEmptyBucketPtr == nullptr) {
                    firstEmptyBucketPtr = &bucket;
                }
            } else if (bucket.fingerprint == fingerprint) {
                if (findedKeyNodePtr != nullptr || bucket.counter <= minValInSS) {
                    bucket.counter++;
                    addFlag = true;
                    toUpdateValue = bucket.counter;
                    break;
                }
                toUpdateValue = max(bucket.counter, toUpdateValue);
            } else {
                if (bucket.counter < counterOfMinBucket) {
                    counterOfMinBucket = bucket.counter;
                    minBucketPtr = &bucket;
                }
            }
        }

        if (addFlag == false) {
            if (firstEmptyBucketPtr != nullptr) {
                firstEmptyBucketPtr->fingerprint = fingerprint;
                firstEmptyBucketPtr->counter = 1;
                toUpdateValue = 1;
            } else if (minBucketPtr != nullptr) {
                double randomVal = (double)rng()/4294967296;
                if (randomVal < pow(DECAY_PARAM, -1.0 * counterOfMinBucket)) {
                    if (counterOfMinBucket == 1) {
                        minBucketPtr->fingerprint = fingerprint;
                        toUpdateValue = 1;
                    } else {
                        minBucketPtr->counter -= 1;
                    }
                }
            }
        }

        if (findedKeyNodePtr != nullptr) {
            if (toUpdateValue > findedKeyNodePtr->parent->counter) {
                ss.SSUpdate(findedKeyNodePtr, toUpdateValue);
            }
        } else {
            if (ss.getSize() < SS_MAX_SIZE or toUpdateValue - minValInSS == 1) {
                ss.SSPush(key, toUpdateValue);
            }
        }
    }


    void getEstimatedFlowSizes(unordered_map<uint8_t *, unsigned int, HashFunc, CmpFunc> &estimatedFlowSizes) {
        SSValueNode *ssValueNode = ss.getFirstValueNode();
        while (ssValueNode != nullptr) {
            uint32_t counter = ssValueNode->counter;
            SSKeyNode *ssKeyNode = ssValueNode->firstKeyNode;

            do {
                uint8_t *tempKey = (uint8_t *) malloc(KEY_SIZE);
                memcpy(tempKey, ssKeyNode->key, KEY_SIZE);
                estimatedFlowSizes[tempKey] = counter;

                ssKeyNode = ssKeyNode->next;
            } while (ssKeyNode != ssValueNode->firstKeyNode);

            ssValueNode = ssValueNode->next;
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
    const uint32_t ssMaxSize =1000;
    const uint32_t colNum = 15712;//200KB
    HeavyKeeper<rowNum,colNum,ssMaxSize> heavyKeeper;

    //output memory info
    double bucketMem = (2 + 18 / 8.0) * rowNum * colNum / 1024.0;
    //suppose there are SS_MAX_SIZE distinct frequencies of top-k flows after processing all items, thus each top-k flow will consume 1 counter, 5 pointers, and a key slot in stream summay's main structure
    //for simplicity, we suppose the hash table consists of SS_MAX_SIZE hash nodes, each of which has 2 pointers (one points to the key slot, the other points to the key node in stream summary)
    //Though we use 32-bit counters, 18 bits is large enough to record the largest flow, thus we use 18 bits as the counter size
    //in 64-bits program, each pointer is 8 bytes, and the key slot is 13 bytes
    double SSMem = ssMaxSize * (KEY_SIZE + 8 * 5 + 8 * 2 + (18 / 8.0)) / 1024.0;
    cout << "bucketMem: " << bucketMem << "KB" << endl;
    cout << "streamSummaryMem: " << SSMem << "KB" << endl;
    cout << "totalMem:" << bucketMem + SSMem << "KB" << endl;
    cout << "*********************" << endl;

    //insert items
    cout << "insert items" << endl;
    clock_t time1 = clock();
    for (unsigned int i = 0; i < actualItemNum; i++) {
        heavyKeeper.insert(keys[i]);
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
    heavyKeeper.getEstimatedFlowSizes(estimatedFlowSizes);

    for (auto iter = estimatedFlowSizes.begin(); iter != estimatedFlowSizes.end(); iter++) {
        estimatedFlowSizesVector.push_back(make_pair(iter->first, iter->second));
    }

    cout << "get top-" << estimatedFlowSizes.size() << " flows" << endl;
    sort(estimatedFlowSizesVector.begin(), estimatedFlowSizesVector.end(), cmpPairFunc);
    unsigned int k = min((unsigned int) ssMaxSize, (unsigned int) estimatedFlowSizes.size());

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
    cout << "Precision: " << (double) TPNumOfTopK / k << endl;

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