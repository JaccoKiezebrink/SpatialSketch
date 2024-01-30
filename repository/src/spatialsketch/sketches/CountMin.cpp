#include "CountMin.h"

#include "MurmurHash3.h"

#include <math.h>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <limits.h>

bool CountMin::isPrime(long prime) {
    if (prime <= 1) {
        return false;
    }
    if (prime <= 3) {
        return true;
    }
    if (prime % 2 == 0 || prime % 3 == 0) {
        return false;
    }
    if (prime % 5 == 0 || prime % 7 == 0) {
        return false;
    }
    if (prime % 11 == 0 || prime % 13 == 0) {
        return false;
    }
    if (prime % 17 == 0 || prime % 19 == 0) {
        return false;
    }
    return true;

}

CountMin::CountMin(float epsilon, float delta, int** hashab) : Sketch() {
    width_ = (int) ceil(std::exp(1) / epsilon);
    repetitions_ = (int) ceil(log(1 / delta));
    counters_ = new int*[repetitions_];

    if (hashab != nullptr) {
        hashab_ = hashab;
        //std::cout << "hashab passed\n";
    }

    for (int i = 0; i < repetitions_; i++) {
        counters_[i] = new int[width_];
        memset(counters_[i], 0, sizeof(int) * width_);  // ensure counters are set to zero

        if (hashab == nullptr) {
            if (i == 0) {
                hashab_ = new int*[repetitions_];
            }
            hashab_[i] = new int[3];
            std::srand((i+1)*7);
            //std:: cout << std::rand() << ", " << std::rand() << ", " << std::rand() << std::endl;
            //hashab_[i][0] = int(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
            //hashab_[i][1] = int(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
            //hashab_[i][0] = int(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
            hashab_[i][0] = std::rand();//long(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
            hashab_[i][1] = std::rand();//long(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
            // generate random large prime number
            long prime = 10000 + ( std::rand() % ( RAND_MAX - 10000 + 1 ) );
            int count = 0;
            while (!isPrime(prime) && count < 1000) {
                prime = 10000 + ( std::rand() % ( RAND_MAX - 10000 + 1 ) );
            }
            
            hashab_[i][2] = prime;
            
            //std::cout << i << ", " <<  hashab[i][0] << " " << hashab[i][1] << std::endl;
        }
    }
}


CountMin::~CountMin() {
    for (int i = 0; i < repetitions_; i++) {
        delete[] counters_[i];
    }
    delete[] counters_;
}


void CountMin::PrintCounters() {
    for (int i = 0; i < repetitions_; i++) {
        for (int j = 0; j < width_; j++) {
            std::cout << counters_[i][j] << " ";
        }
        std::cout << std::endl;
    }
}


/*void CountMin::SetupHashes(int** hashab) {
    hashab = new int*[depth_];

    for (int i = 0; i < depth_; i++) {
        hashab[i] = new int[2];

        std::srand((i+1)*7);
        //std:: cout << std::rand() << ", " << std::rand() << ", " << std::rand() << std::endl;
        hashab[i][0] = int(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
        hashab[i][1] = int(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
        //std::cout << i << ", " <<  hashab[i][0] << " " << hashab[i][1] << std::endl;
    }
}*/

// Given an item id, compute the d hashes for it to be reused
inline int CountMin::GetItemHashes(long id, uint* hashes) {
    //memset(hashes, 0, sizeof(uint) * depth_);
    for (int i = 0; i < repetitions_; i++) {
        /*uint hash;
        MurmurHash3_x86_32(&id, sizeof(int), i, &hash);
        hashes[i] = hash % width_;*/
        //hashes[i] = (hashab_[i][0] * id + hashab_[i][1]) % width_;
        hashes[i] = ((hashab_[i][0] * id + hashab_[i][1]) % hashab_[i][2] % width_ + width_) % width_;
    }
    return repetitions_;
}


// Insert item in count-min, hashes not given
inline void CountMin::Insert(long id) {
    for (int i = 0; i < repetitions_; i++) {
        /*uint hash;
        MurmurHash3_x86_32(&id, sizeof(int), i, &hash);
        counters_[i][hash % width_] += count;*/
        //counters_[i][(hashab_[i][0] * id + hashab_[i][1]) % width_] += 1;
        counters_[i][((hashab_[i][0] * id + hashab_[i][1]) % hashab_[i][2] % width_ + width_) % width_] += 1;
    }
}


// Insert item in count-min, hashes already given
inline void CountMin::Insert(long id, uint* hashes) {
    for (int i = 0; i < repetitions_; i++) {
        counters_[i][hashes[i]] += 1;
    }
}

// Insert item in count-min, hashes not given

inline void CountMin::Insert(long id, int count) {
    for (int i = 0; i < repetitions_; i++) {
        /*uint hash;
        MurmurHash3_x86_32(&id, sizeof(int), i, &hash);
        counters_[i][hash % width_] += count;*/
        //counters_[i][(hashab_[i][0] * id + hashab_[i][1]) % width_] += count;
        counters_[i][((hashab_[i][0] * id + hashab_[i][1]) % hashab_[i][2] % width_ + width_) % width_] += count;
    }
}


// Insert item in count-min, hashes already given
inline void CountMin::Insert(long id, int count, uint* hashes) {
    for (int i = 0; i < repetitions_; i++) {
        counters_[i][hashes[i]] += count;
    }
}

// Query the respective d counters for the given item id, hashes not given
std::vector<int> CountMin::QueryCounters(long id) {
    std::vector<int> counters = std::vector<int>(repetitions_);
    for (int i = 0; i < repetitions_; i++) {
        /*uint hash;
        MurmurHash3_x86_32(&id, sizeof(int), i, &hash);
        counters[i] = counters_[i][hash % width_];*/
        //counters[i] = counters_[i][(hashab_[i][0] * id + hashab_[i][1]) % width_]; //[hash % width_];
        counters[i] = counters_[i][((hashab_[i][0] * id + hashab_[i][1]) % hashab_[i][2] % width_ + width_) % width_];
    }

    return counters;
}

// Query the respective d counters for the given item id, hashes already given
std::vector<int> CountMin::QueryCounters(long id, uint* hashes) {
    std::vector<int> counters = std::vector<int>(repetitions_);
    for (int i = 0; i < repetitions_; i++) {
        counters[i] = counters_[i][hashes[i]];
    }

    return counters;
}

int CountMin::QueryItem(long id) {
    std::vector<int> counters;
    counters = QueryCounters(id);
    return *std::min_element(counters.begin(), counters.end());
}

int CountMin::QueryItem(long id, uint* hashes) {
    std::vector<int> counters;
    counters = QueryCounters(id, hashes);
    return *std::min_element(counters.begin(), counters.end());
}

void CountMin::Merge(Sketch *cm) {
    for (int i = 0; i < repetitions_; i++) {
        for (int j = 0; j < width_; j++) {
            counters_[i][j] += cm->counters_[i][j];
        }
    }
}


long CountMin::L2Estimate() {
    long min = LONG_MAX;
    for (int i = 0; i < repetitions_; i++) {
        long sum = 0;
        for (int j = 0; j < width_; j++) {
            sum += (long) std::pow(counters_[i][j], 2);
        }
        min = std::min(min, sum);
    }
    return min;
}
