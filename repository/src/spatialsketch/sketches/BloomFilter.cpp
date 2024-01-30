#include "BloomFilter.h"
#include <math.h>
#include <iostream>

bool BloomFilter::isPrime(long prime) {
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

BloomFilter::BloomFilter(float delta,  int N, long** hashab) : Sketch() {
    bits_ = (long) ceil(-N * log(delta) / (log(2) * log(2)));
    
    repetitions_ = (int) ceil(- log(delta) / log(2));

    // std::cout << "BF bits: " << bits_ << std::endl;
    // std::cout << "BF repetitions: " << repetitions_ << std::endl;


    bloomFilter = std::vector<bool>(bits_, false);

    if (hashab != nullptr) {
        hashab_ = hashab;
        //std::cout << "hashab passed\n";
    }

    for (int i = 0; i < repetitions_; i++) {
        if (hashab == nullptr) {
            if (i == 0) {
                hashab_ = new long*[repetitions_];
            }
            hashab_[i] = new long[3];
            std::srand((i+1)*7);
            //std:: cout << std::rand() << ", " << std::rand() << ", " << std::rand() << std::endl;
            // hashab_[i][0] = int(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
            // hashab_[i][1] = int(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
            
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

int BloomFilter::GetItemHashes(long id, uint* hashes) {
    for (int i = 0; i < repetitions_; i++) {
        /*uint hash;
        MurmurHash3_x86_32(&id, sizeof(int), i, &hash);
        hashes[i] = hash % width_;*/
        hashes[i] = (int) (((hashab_[i][0] * id + hashab_[i][1]) % hashab_[i][2] % bits_) + bits_) % bits_;
    }
    return repetitions_;
}

void BloomFilter::Insert(long id) {   
    //std::cout << "BF: Inserting " << id << " with repetitions: " << repetitions_ << std::endl;
    int first_index = (int) (((hashab_[0][0] * id + hashab_[0][1]) % hashab_[0][2] % bits_) + bits_) % bits_;
    for (int i = 0; i < repetitions_; i++) {
        //int index = sha.hash(id, i); % bits_;
        int index = (int) (((hashab_[i][0] * id + hashab_[i][1]) % hashab_[i][2] % bits_) + bits_) % bits_;
        // if (index < 0) {
        //     std::cout << "index < 0: " << index << std::endl;
        // } else {
        //     std::cout << "index: " << index << std::endl;
        // }
        //bool before = bloomFilter[index];
        //std::cout << "BF: Inserting " << id << " at index " << index << " previously " << before << std::endl;
        bloomFilter[index] = true;
        //bool after = bloomFilter[index];
        // if (before && !after) {
        //     throw "BF: Insertion error";
        // } else if (!before && after) {
        //     std::cout << "BF: Insertion success" << std::endl;
        // } else if (!before && !after) {
        //     throw "BF: Insertion error";
        // } else if (before && after) {
        //     std::cout << "BF: Insertion success" << std::endl;
        // }
    }

    // check if first index is still true
    if (!bloomFilter[first_index]) {
        throw "BF: Insertion error";
    }
}

void BloomFilter::Insert(long id, uint* hashes) {    
    for (int i = 0; i < repetitions_; i++) {
        bloomFilter[hashes[i]] = true;

    }
}

void BloomFilter::Insert(long id, int count) {    
    Insert(id);
}

void BloomFilter::Insert(long id, int count, uint* hashes) {    
    Insert(id, hashes);
}



int BloomFilter::QueryItem(long id) {
    for (int i = 0; i < repetitions_; i++) {
        int index =(int) (((hashab_[i][0] * id + hashab_[i][1]) % hashab_[i][2] % bits_) + bits_ ) % bits_;
        if (!bloomFilter[index]) {
            return 0;
        }
    }
    return 1;
}

int BloomFilter::QueryItem(long id, uint* hashes) {
    throw "BF does not support hashes";
    //return QueryItem(id);
}
