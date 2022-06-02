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

//the key size cannot be simply changed since other codes are designed for 13-byte IDs
#define KEY_SIZE 13
#define MASK_26BITS 0x3FFFFFF

//the exp part and coe part share 16 bits, and bot of them cannot be less than 1.
#define EXP_PART_BITS 2
#define COE_PART_BITS 14
//these are corresponding definitions to reduce calculation
#define EXP_PART_MASK 3
#define COE_PART_DEC_BITS 13
#define COE_PART_DEC_MASK 0x1FFF
#define COE_PART_BASEVAL 0x2000

#define EXTRA_BITS_NUM 2
//a pair of randomized generated modular inverses
#define MI_A 1635265620366789
#define MI_A_INV 4417919481176333
// MI_MASK is pow(2,52)-1 here, used for replacing the mod operation
#define MI_MASK 4503599627370495

static mt19937 rng(time(0));

bool cmpPairFunc(pair<uint8_t*, int>p1, pair<uint8_t*, int>p2)
{
    return p1.second > p2.second;
}

//for recording information in unordered_map
struct CmpFunc {
    bool operator()(const uint8_t* keyA, const uint8_t* keyB) const {
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

template<uint32_t BUCKET_NUM,uint32_t LEFT_PART_BITS,uint32_t CELL_NUM_H,uint32_t CELL_NUM_L>
class SKETCH{
private:
    struct node {uint16_t FP=0; uint16_t C=0;} bucketArray[BUCKET_NUM][CELL_NUM_H + CELL_NUM_L];
    uint64_t* auxiliaryList;
public:
    SKETCH(){
        uint32_t auxiliaryListWordNum = int(ceil(BUCKET_NUM * CELL_NUM_H * (LEFT_PART_BITS+EXTRA_BITS_NUM) / 64.0));
        auxiliaryList = (uint64_t*)calloc(auxiliaryListWordNum, 8);
    }

    void getEstimatedFlowSizes(unordered_map<uint8_t*, unsigned int, HashFunc, CmpFunc>& estimatedFlowSizes) {
        for (uint32_t bucketIdx = 0; bucketIdx < BUCKET_NUM; bucketIdx++) {
            //traverse heavy cells
            for (int i = 0; i < CELL_NUM_H; i++) {
                uint16_t fp = bucketArray[bucketIdx][i].FP;
                uint16_t counterOfCurrentCell = bucketArray[bucketIdx][i].C;
                uint32_t counterRealVal = getRealValOfTwoMoldActiveCounter(counterOfCurrentCell);

                uint64_t leftPart[2] = { 0 };
                uint32_t slotIndex = bucketIdx * CELL_NUM_H + i;
                getLeftPart(slotIndex, leftPart);

                uint8_t* reverseKey = (uint8_t*)malloc(KEY_SIZE);
                combineKey(reverseKey, bucketIdx, fp, leftPart);

                estimatedFlowSizes[reverseKey] = counterRealVal;
            }
        }
    }

    void insert(uint8_t* key) {
        uint32_t bucketIdx = 0;
        uint16_t fp = 0;
        uint64_t leftPart[2] = { 0 };
        divideKey(key, bucketIdx, fp, leftPart);

        int matchedCellIdx = -1;
        uint32_t matchedCellCounter = -1;

        int smallestHeavyCellIdx = -1;
        uint32_t smallestHeavyCellCounter = -1;
        uint16_t smallestHeavyCellFP = 0;

        int smallestCellIdx = -1;
        uint32_t smallestCellCounter = -1;

        uint16_t fingerprintOfCurrentCell;
        uint16_t counterOfCurrentCell;

        //traverse heavy cells
        for (int i = 0; i < CELL_NUM_H; i++) {
            fingerprintOfCurrentCell = bucketArray[bucketIdx][i].FP;
            counterOfCurrentCell = bucketArray[bucketIdx][i].C;
            if(counterOfCurrentCell==0){
                bucketArray[bucketIdx][i].FP=fp;
                bucketArray[bucketIdx][i].C=1;

                uint32_t  slotIndex = bucketIdx * CELL_NUM_H + i;
                setLeftPart(slotIndex, leftPart);
                return;
            }

            if (fp == fingerprintOfCurrentCell && counterOfCurrentCell > 0) {
                matchedCellIdx = i;
                matchedCellCounter = counterOfCurrentCell;
                break;
            }

            if (counterOfCurrentCell < smallestHeavyCellCounter) {
                smallestHeavyCellIdx = i;
                smallestHeavyCellFP = fingerprintOfCurrentCell;
                smallestHeavyCellCounter = counterOfCurrentCell;
            }

            if (counterOfCurrentCell < smallestCellCounter) {
                smallestCellIdx = i;
                smallestCellCounter = counterOfCurrentCell;
            }
        }

        //traverse light cells
        if (matchedCellIdx < 0) {
            for (int i = CELL_NUM_H; i < CELL_NUM_H + CELL_NUM_L; i++) {
                fingerprintOfCurrentCell = bucketArray[bucketIdx][i].FP;
                counterOfCurrentCell = bucketArray[bucketIdx][i].C;

                if(counterOfCurrentCell==0){
                    bucketArray[bucketIdx][i].FP=fp;
                    bucketArray[bucketIdx][i].C=1;
                    return;
                }
                if (fp == fingerprintOfCurrentCell && counterOfCurrentCell > 0) {
                    matchedCellIdx = i;
                    matchedCellCounter = counterOfCurrentCell;
                    break;
                }

                if (counterOfCurrentCell < smallestCellCounter) {
                    smallestCellIdx = i;
                    smallestCellCounter = counterOfCurrentCell;
                }
            }
        }

        if (matchedCellIdx < 0) {//do not find a matched cell that has the same fingerprint
            uint32_t smallestCellCounterVal = getRealValOfTwoMoldActiveCounter(smallestCellCounter);

            //replace under probabiltiy
            double prob = 1.0 / smallestCellCounterVal;
            double randomVal = (double)rng()/4294967296;
            if (randomVal < prob) {
                bucketArray[bucketIdx][smallestCellIdx].FP=fp;
                if (smallestCellIdx < CELL_NUM_H) {
                    uint32_t  slotIndex = bucketIdx * CELL_NUM_H + smallestCellIdx;
                    setLeftPart(slotIndex, leftPart);
                }
            }
        }
        else {
            uint16_t newCounter = incTwoMoldActiveCounter(matchedCellCounter);
            if (newCounter != matchedCellCounter) {
                if (matchedCellIdx >= CELL_NUM_H && newCounter > smallestHeavyCellCounter) {
                    bucketArray[bucketIdx][matchedCellIdx].FP = smallestHeavyCellFP;
                    bucketArray[bucketIdx][matchedCellIdx].C = smallestHeavyCellCounter;
                    bucketArray[bucketIdx][smallestHeavyCellIdx].FP = fp;
                    bucketArray[bucketIdx][smallestHeavyCellIdx].C = newCounter;

                    uint32_t slotIndex = bucketIdx * CELL_NUM_H + smallestHeavyCellIdx;
                    setLeftPart(slotIndex, leftPart);
                    return;
                } else {
                    bucketArray[bucketIdx][matchedCellIdx].C = newCounter;
                }
            }

            if(matchedCellIdx < CELL_NUM_H && newCounter>=100){
                uint32_t newCounterVal= getRealValOfTwoMoldActiveCounter(newCounter);

                double prob = 1.0 / newCounterVal;
                double randomVal = (double)rng()/4294967296;
                if(randomVal<prob){
                    uint32_t  slotIndex = bucketIdx * CELL_NUM_H + matchedCellIdx;
                    uint64_t targetleftPart[2] = { 0 };
                    uint8_t ExtraCounter=getLeftPart(slotIndex, targetleftPart);
                    if(memcmp(leftPart,targetleftPart,16)!=0){
                        if(ExtraCounter>0){
                            setLeftPartCounter(slotIndex, ExtraCounter-1);
                        }else{
                            setLeftPart(slotIndex, leftPart);
                        }
                    }else{
                        if(ExtraCounter!=(1<<(EXTRA_BITS_NUM))-1){
                            setLeftPartCounter(slotIndex, ExtraCounter+1);
                        }
                    }
                }
            }
        }
    }

    ~SKETCH(){
        free(auxiliaryList);
    }

private:
    uint8_t getLeftPart(uint32_t slotIndex, uint64_t* leftPart) {
        uint8_t counter=0;

        unsigned int bitIdx = slotIndex * (LEFT_PART_BITS+EXTRA_BITS_NUM);
        unsigned int slotWordIdx = bitIdx / 64;
        unsigned int slotBitIdxInWord = bitIdx % 64;

        unsigned int extractedBitsNum = 0;
        unsigned int LPWordIdx = 0;
        unsigned int LPBitInWord = 0;

        while (extractedBitsNum < LEFT_PART_BITS) {
            unsigned int toExtractBitsNum = min(LEFT_PART_BITS - extractedBitsNum, 64 - LPBitInWord);
            toExtractBitsNum = min(toExtractBitsNum, 64 - slotBitIdxInWord);

            uint64_t extractPart = 0;
            uint64_t extractPartMask = 0;
            if (toExtractBitsNum == 64) {
                extractPart = auxiliaryList[slotWordIdx];
            }
            else {
                extractPartMask = ((uint64_t)1 << toExtractBitsNum) - 1;
                extractPart = (auxiliaryList[slotWordIdx] >> slotBitIdxInWord) & extractPartMask;
            }
            if (LPBitInWord == 0) {
                leftPart[LPWordIdx] = 0;
            }
            leftPart[LPWordIdx] += extractPart << LPBitInWord;

            bitIdx += toExtractBitsNum;
            slotWordIdx = bitIdx / 64;
            slotBitIdxInWord = bitIdx % 64;

            extractedBitsNum += toExtractBitsNum;
            LPWordIdx = extractedBitsNum / 64;
            LPBitInWord = extractedBitsNum % 64;
        }

        extractedBitsNum = 0;
        while (extractedBitsNum < EXTRA_BITS_NUM) {
            unsigned int toExtractBitsNum = min(EXTRA_BITS_NUM - extractedBitsNum, 64 - slotBitIdxInWord);

            uint64_t extractPart = 0;
            uint64_t extractPartMask = 0;

            extractPartMask = ((uint64_t)1 << toExtractBitsNum) - 1;
            extractPart = (auxiliaryList[slotWordIdx]>> slotBitIdxInWord) & extractPartMask;

            counter += extractPart << extractedBitsNum;

            bitIdx += toExtractBitsNum;
            slotWordIdx = bitIdx / 64;
            slotBitIdxInWord = bitIdx % 64;

            extractedBitsNum += toExtractBitsNum;
        }
        return counter;
    }

    void setLeftPart(uint32_t slotIndex, uint64_t* leftPart) {
        unsigned int bitIdx = slotIndex * (LEFT_PART_BITS+EXTRA_BITS_NUM);
        unsigned int slotWordIdx = bitIdx / 64;
        unsigned int slotBitIdxInWord = bitIdx % 64;

        unsigned int extractedBitsNum = 0;
        unsigned int LPWordIdx = 0;
        unsigned int LPBitInWord = 0;

        while (extractedBitsNum < LEFT_PART_BITS) {
            unsigned int toExtractBitsNum = min(LEFT_PART_BITS - extractedBitsNum, 64 - LPBitInWord);
            toExtractBitsNum = min(toExtractBitsNum, 64 - slotBitIdxInWord);

            uint64_t extractPart = 0;
            uint64_t extractPartMask = 0;
            if (toExtractBitsNum == 64) {
                extractPart = leftPart[LPWordIdx];

                auxiliaryList[slotWordIdx] = extractPart;
            }
            else {
                extractPartMask = ((uint64_t)1 << toExtractBitsNum) - 1;
                extractPart = (leftPart[LPWordIdx] >> LPBitInWord) & extractPartMask;

                auxiliaryList[slotWordIdx] &= ~(extractPartMask << slotBitIdxInWord);//clear the corresponding bits
                auxiliaryList[slotWordIdx] += extractPart << slotBitIdxInWord;
            }

            bitIdx += toExtractBitsNum;
            slotWordIdx = bitIdx / 64;
            slotBitIdxInWord = bitIdx % 64;

            extractedBitsNum += toExtractBitsNum;
            LPWordIdx = extractedBitsNum / 64;
            LPBitInWord = extractedBitsNum % 64;
        }
    }
    void setLeftPartCounter(uint32_t slotIndex,uint8_t counter) {
        unsigned int bitIdx = slotIndex * (LEFT_PART_BITS+EXTRA_BITS_NUM)+LEFT_PART_BITS;

        unsigned int slotWordIdx = bitIdx / 64;
        unsigned int slotBitIdxInWord = bitIdx % 64;

        unsigned int extractedBitsNum = 0;

        while (extractedBitsNum < EXTRA_BITS_NUM) {
            unsigned int toExtractBitsNum = min(EXTRA_BITS_NUM - extractedBitsNum, 64 - slotBitIdxInWord);

            uint64_t extractPart = 0;
            uint64_t extractPartMask = 0;

            extractPartMask = ((uint64_t)1 << toExtractBitsNum) - 1;
            extractPart = (counter>>extractedBitsNum) & extractPartMask;

            auxiliaryList[slotWordIdx] &= ~(extractPartMask << slotBitIdxInWord);//clear the corresponding bits
            auxiliaryList[slotWordIdx] += extractPart << slotBitIdxInWord;


            bitIdx += toExtractBitsNum;
            slotWordIdx = bitIdx / 64;
            slotBitIdxInWord = bitIdx % 64;

            extractedBitsNum += toExtractBitsNum;
        }
    }

    //get the real value represented by the two mode active counter
    uint32_t getRealValOfTwoMoldActiveCounter(uint16_t counter) {
        if (counter > 0x8000) {//in active mode
            uint32_t coePart = (counter & COE_PART_DEC_MASK) + COE_PART_BASEVAL;
            uint32_t expPart = (counter >> COE_PART_DEC_BITS) & EXP_PART_MASK;

            uint32_t value = coePart << (EXP_PART_BITS + expPart);
            return value;
        }
        else {//in normal mode
            return counter;
        }
    }

    uint16_t incTwoMoldActiveCounter(uint16_t counter) {
        if (counter >= 0x8000) {//in active mode
            uint32_t coePart = (counter & COE_PART_DEC_MASK);
            uint32_t expPart = (counter >> COE_PART_DEC_BITS) & EXP_PART_MASK;

            double randomVal = (double)rng()/4294967296;
            if (randomVal * ((uint32_t)1 << (EXP_PART_BITS + expPart)) < 1.0) {
                if (coePart < COE_PART_DEC_MASK) {//the coe part can be directly increased by one without overflow
                    coePart += 1;
                }
                else {
                    if (expPart < EXP_PART_MASK) {
                        coePart = 0;
                        expPart += 1;
                    }
                }

                uint16_t newCounter = 0x8000 + (expPart << COE_PART_DEC_BITS) + coePart;
                return newCounter;
            }
            else {
                return counter;
            }
        }
        else {//in normal mode
            return counter + 1;
        }
    }

    //designed for five-tuple key
    void divideKey(uint8_t key[KEY_SIZE], uint32_t& index, uint16_t& fingerprint, uint64_t leftPart[2]) {
        memcpy(leftPart, key, KEY_SIZE);

        //encode the key using modular inverse
        uint64_t part1 = leftPart[0] & MI_MASK;
        uint64_t part2 = (leftPart[1] << 12) + (leftPart[0] >> 52);
        part1 = (part1 * MI_A) & MI_MASK;
        part2 = (part2 * MI_A) & MI_MASK;

        uint32_t tempParts[2] = { 0 };
        tempParts[0] = (uint32_t)(part1 & MASK_26BITS) ^ (uint32_t)(part1 >> 26) ^ (uint32_t)(part2 & MASK_26BITS);
        tempParts[1] = tempParts[0] ^ (uint32_t)(part2 >> 26);
        tempParts[1] ^= tempParts[1] >> 13;

        leftPart[0] = part1 + ((uint64_t)(tempParts[0] & 0xFFF) << 52);
        leftPart[1] = (tempParts[1]) + (((uint64_t)tempParts[0] & (~0xFFF)) <<14);

        index = leftPart[1] % BUCKET_NUM;
        leftPart[1] /= BUCKET_NUM;
        fingerprint = leftPart[1] & (0xffff);
        leftPart[1] >>= 16;
    }

    //designed for five-tuple key
    void combineKey(uint8_t key[KEY_SIZE], uint32_t index, uint16_t fingerprint, uint64_t leftPart[2]) {
        leftPart[1] = ((leftPart[1] << 16) + fingerprint) * BUCKET_NUM + index;
        uint32_t tempParts[2] = { 0 };
        tempParts[0] = (leftPart[0] >> 52) + ((leftPart[1] >> 26) << 12);
        tempParts[1] = leftPart[1] & 0x3FFFFFF;

        tempParts[1] = (tempParts[1] & 0x1FFF) ^ (tempParts[1] >> 13) + (tempParts[1] & (~0x1FFF));

        uint64_t part1 = leftPart[0] & MI_MASK;
        uint64_t part2 = 0;

        part2 += tempParts[1] ^ tempParts[0];
        part2<<= 26;
        part2 += tempParts[0] ^ (uint32_t)(part1 & MASK_26BITS) ^ (uint32_t)(part1 >> 26);

        part1 = (part1 * MI_A_INV) & MI_MASK;
        part2 = (part2 * MI_A_INV) & MI_MASK;

        leftPart[0] = part1 + ((part2 & 0xfff) << 52);
        leftPart[1] = part2 >> 12;
        memcpy(key, leftPart, 13);
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
    cout << "prepare algorithm "<< endl;
    //parameters for 200KB memory on finding top-1000 flows
    const uint32_t bucketNum=3324;
    const uint32_t leftPartBits=77;//this value have to equal (KEY_SIZE*8)-16-floor(log2(bucketNum));
    const uint32_t cellNum_H=3;
    const uint32_t cellNum_L=5;

    uint32_t tempK=1000;
    SKETCH<bucketNum,leftPartBits,cellNum_H,cellNum_L> sketch;

    //output memory info
    uint32_t temp=KEY_SIZE*8-16-floor(log2(bucketNum));
    double bucketMem = bucketNum * (cellNum_H + cellNum_L) * 4 / 1024.0;
    uint32_t auxiliaryListWordNum = int(ceil(bucketNum * cellNum_H  *(leftPartBits+EXTRA_BITS_NUM) / 64.0));
    double auxiliaryListMem = auxiliaryListWordNum * 8 / 1024.0;    cout << "bucketMem: " << bucketMem << "KB" << endl;
    cout << "auxiliaryListMem: " << auxiliaryListMem << "KB" << endl;
    cout << "totalMem:" << bucketMem + auxiliaryListMem << "KB" << endl;
    cout << "*********************" << endl;

    //insert items
    cout << "insert items" << endl;
    clock_t time1 = clock();
    for (unsigned int i = 0; i < actualItemNum; i++) {
        sketch.insert(keys[i]);
    }
    clock_t time2 = clock();

    //calculate throughput
    double numOfSeconds = (double)(time2 - time1) / CLOCKS_PER_SEC;//the seconds using to insert items
    double throughput = (actualItemNum / 1000000.0) / numOfSeconds;
    cout << "use " << numOfSeconds << " seconds" << endl;
    cout << "throughput: " << throughput << " Mpps" << ", each insert operation uses " << 1000.0 / throughput << " ns" << endl;
    cout << "*********************" << endl;

    //calculate precision, ARE, AAE

    //get sorted actual flow sizes and sorted estimated flow sizes
    vector<pair<uint8_t*, unsigned int>> actualFlowSizesVector;
    for (auto iter = actualFlowSizes.begin(); iter != actualFlowSizes.end(); iter++) {
        actualFlowSizesVector.push_back(make_pair(iter->first, iter->second));
    }
    sort(actualFlowSizesVector.begin(), actualFlowSizesVector.end(), cmpPairFunc);
    vector<pair<uint8_t*, unsigned int>> estimatedFlowSizesVector;
    unordered_map<uint8_t*, unsigned int, HashFunc, CmpFunc> estimatedFlowSizes;
    sketch.getEstimatedFlowSizes(estimatedFlowSizes);

    for (auto iter = estimatedFlowSizes.begin(); iter != estimatedFlowSizes.end(); iter++) {
        estimatedFlowSizesVector.push_back(make_pair(iter->first, iter->second));
    }

    sort(estimatedFlowSizesVector.begin(), estimatedFlowSizesVector.end(), cmpPairFunc);
    cout << "actually recovered " << estimatedFlowSizes.size() << " flows" << endl;
    unsigned int k = min((unsigned int)tempK, (unsigned int)estimatedFlowSizes.size());


    //get acutal top-k flows
    unordered_map<uint8_t*, unsigned int, HashFunc, CmpFunc> actualTopKFlowSizes;
    for (unsigned int i = 0; i < k; i++) {
        unsigned int actualSize = actualFlowSizesVector[i].second;
        actualTopKFlowSizes[actualFlowSizesVector[i].first] = actualSize;
    }
    cout << "Top-" << k << " flow detection Metrics" << endl;

    //top-k precision
    unsigned int TPNumOfTopK = 0;//True positive num
    for (int i = 0; i < k; i++) {
        uint8_t* key = estimatedFlowSizesVector[i].first;
        if (actualTopKFlowSizes.find(key) != actualTopKFlowSizes.end()) {
            TPNumOfTopK++;
        }
    }
    cout << "Precision: " << (double)TPNumOfTopK / k << endl;

    //ARE,AAE of the reported top-k flows
    double totalREOfTopK = 0;//total relative error of top-k flows
    double totalAEOfTopK = 0;
    for (int i = 0; i < k; i++) {
        uint8_t* key = estimatedFlowSizesVector[i].first;
        unsigned int estimatedSize = estimatedFlowSizesVector[i].second;
        unsigned int actualSize = 0;
        if (actualFlowSizes.find(key) != actualFlowSizes.end()) {
            actualSize = actualFlowSizes[key];
        }
        totalAEOfTopK += abs((double)estimatedSize - actualSize);
        double rE = abs((double)estimatedSize - actualSize) / actualSize;
        totalREOfTopK += rE;
    }
    cout << "ARE of reported top-k flows: " << totalREOfTopK / k << endl;
    cout << "AAE of reported top-k flows: " << totalAEOfTopK / k << endl;

    //release resources
    for (auto iter = estimatedFlowSizes.begin(); iter != estimatedFlowSizes.end(); iter++) {
        free(iter->first);
    }
    for(unsigned int i=0;i<actualItemNum;i++){
        free(keys[i]);
    }
    free(keys);
    return 0;
}