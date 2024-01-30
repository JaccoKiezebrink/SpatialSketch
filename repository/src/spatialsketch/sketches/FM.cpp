#include "FM.h"
#include <iostream>
#include <math.h>
//#include "Utils.h" 
bool FM::isPrime(long prime) {
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
FM::FM(float epsilon, float delta, long** hashab) : Sketch() {
  
    repetitions_ = (int) ceil( 1/(epsilon* epsilon) * log(1/delta));
    //how does delta influence parameters? b_ is already the domain size. 
    // m bits in register in structure where its easy to get index of first 1;
    r_bitmap = std::vector<std::vector<bool> >(repetitions_, std::vector<bool>(b_, false));

    if (hashab != nullptr) {
        hashab_ = hashab;
        //std::cout << "hashab passed" << hashab_[0][0] << " " << hashab_[0][1] << std::endl}
    }

    for (int i = 0; i < repetitions_; i++) {
        if (hashab == nullptr) {
            if (i == 0) {
                hashab_ = new long*[repetitions_];
            }
            hashab_[i] = new long[3];
            std::srand((i+1)*7);
            //std:: cout << std::rand() << ", " << std::rand() << ", " << std::rand() << std::endl;
            hashab_[i][0] = std::rand();//long(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
            hashab_[i][1] = std::rand();//long(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
            // generate random large prime number
            long prime = 10000 + ( std::rand() % ( RAND_MAX - 10000 + 1 ) );
            int count = 0;
            while (!isPrime(prime) && count < 1000) {
                prime = 10000 + ( std::rand() % ( RAND_MAX - 10000 + 1 ) );
            }
            
            hashab_[i][2] = prime; //long(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
        }
    }

}
FM::~FM() {}
int Query(dyadic1D dyad) {
    throw "FM does not support dyadic query";
}
void FM::PrintHashab() {
    long initHashab0 = hashab_[0][0];
    long initHashab1 = hashab_[0][1];
    for (int i = 1; i < repetitions_; i++) {
        if (i % 50 == 0) {
            std::cout << " same hash function for " << i << " repetitions out of " << repetitions_ << std::endl;
        }
        if (hashab_[i][0] != initHashab0 || hashab_[i][1] != initHashab1) {
            std::cout << "hashab not the same, namely: " << hashab_[i][0] << " " << hashab_[i][1] << " instead of " << initHashab0 << " " << initHashab1 << " at repetition " << i << std::endl;
            break;
        }

    }
}

inline int trailingZeros(int x) {
    if (x==0) {
        return 1;
    }
    int count = 0;
    while ((x & 1) == 0) {
        count++;
        x >>= 1;
    }
    return count;
}

inline int FM::GetItemHashes(long id, long* hashes) {
    int index;
    for (int i = 0; i < repetitions_; i++) {
        // compute hash
        index = trailingZeros(((hashab_[i][0] * id + hashab_[i][1])% hashab_[i][2] % M_ + M_) % M_);
        hashes[i] = index;
    }
    return repetitions_;
}

void FM::Insert(long id) {
    // compute one hash x = h(v), where we hash int to binary domain
    int index;
    new_insert_ = false;
    for (int i = 0; i < repetitions_; i++) {
        // compute hash
        index = trailingZeros(((hashab_[i][0] * id + hashab_[i][1])% hashab_[i][2] % M_ + M_) % M_);
        if (r_bitmap[i][index] == false) {
            r_bitmap[i][index] = true;
            new_insert_ = true;
        }
    }
}

void FM::Insert(long id, int count) {
    Insert(id);
}

inline void FM::Insert(long id, long* hashes) {
    new_insert_ = false;
    for (int i = 0; i < repetitions_; i++) {
        // compute hash        
        if (!r_bitmap[i][hashes[i]]) {
            r_bitmap[i][hashes[i]] = true;
            new_insert_ = true;
        }
    }
}

inline void FM::Insert(long id, int count, uint* hashes) {
    Insert(id, hashes);
}

int FM::QueryItem(long id) {
    int sum = 0;
    for (int i = 0; i < repetitions_; i++) {
        int leftMostZero;
        for (int j = 0; j < b_; j++) {
            if (r_bitmap[i][j] == false) {
                leftMostZero = j;
                break;
            }
        }
        sum += leftMostZero;
    }
    return (int) pow(2, (float) sum / repetitions_) * 1.2928;
}

int FM::QueryItem(long id, uint* hashes) {
    return FM::QueryItem(id);
}

void MergeFM(FM &a, Sketch *b) {
    // first check if hashab is the same
    // bool sameHashab = true;
    // for (int i = 0; i < a.repetitions_; i++) {
    //     if (a.hashab_[i][0] != b->hashab_[i][0] || a.hashab_[i][1] != b->hashab_[i][1]) {
    //         sameHashab = false;
    //         break;
    //     }
    // }
    // if (!sameHashab) {
    //     throw "Hashab not the same, FMs cannot be merged.";
    // }
    
    for (int i = 0; i < a.repetitions_; i++) {
        for (int j = 0; j < a.b_; j++) {
            a.r_bitmap[i][j] = a.r_bitmap[i][j] || b->r_bitmap[i][j];
        }
    }
}
// FM MergeFM(FM &a, Sketch &b) {
//     FM merged = FM();
//     for (int i = 0; i < a.repetitions_; i++) {
//         for (int j = 0; j < a.b_; j++) {
//             merged.r_bitmap[i][j] = a.r_bitmap[i][j] || b.r_bitmap[i][j];
//         }
//     }
//     return merged;
// }

// #include "FM.h"
// #include <iostream>
// #include <math.h>

// FM::FM(float epsilon, float delta, long** hashab) : Sketch() {
  
//     repetitions_ = (int) ceil( 1/(epsilon* epsilon) * log(1/delta));
//     //how does delta influence parameters? b_ is already the domain size. 
//     // m bits in register in structure where its easy to get index of first 1;
//     r_bitmap = std::vector<std::vector<bool> >(repetitions_, std::vector<bool>(b_, false));

//     if (hashab != nullptr) {
//         hashab_ = hashab;
//         //std::cout << "hashab passed" << hashab_[0][0] << " " << hashab_[0][1] << std::endl;
//     }

//     for (int i = 0; i < repetitions_; i++) {
//         if (hashab == nullptr) {
//             if (i == 0) {
//                 hashab_ = new long*[repetitions_];
//             }
//             hashab_[i] = new long[2];
//             std::srand((i+1)*7);
//             //std:: cout << std::rand() << ", " << std::rand() << ", " << std::rand() << std::endl;
//             hashab_[i][0] = long(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
//             hashab_[i][1] = long(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
//             //std::cout << i << ", " <<  hashab[i][0] << " " << hashab[i][1] << std::endl;
//         }
//     }

// }
// FM::~FM() {}

// void FM::PrintHashab() {
//     long initHashab0 = hashab_[0][0];
//     long initHashab1 = hashab_[0][1];
//     for (int i = 1; i < repetitions_; i++) {
//         if (i % 50 == 0) {
//             std::cout << " same hash function for " << i << " repetitions out of " << repetitions_ << std::endl;
//         }
//         if (hashab_[i][0] != initHashab0 || hashab_[i][1] != initHashab1) {
//             std::cout << "hashab not the same, namely: " << hashab_[i][0] << " " << hashab_[i][1] << " instead of " << initHashab0 << " " << initHashab1 << " at repetition " << i << std::endl;
//             break;
//         }

//     }
// }

// int trailingZeros(int x) {
//     if (x==0) {
//         return 1;
//     }
//     int count = 0;
//     while ((x & 1) == 0) {
//         count++;
//         x >>= 1;
//     }
//     return count;
// }

// int FM::GetItemHashes(long id, uint* hashes) {
//     for (int i = 0; i < repetitions_; i++) {
//         // compute hash
//         int h = (int) ((hashab_[i][0] * id + hashab_[i][1]) % M_ + M_) % M_;
//         int index = trailingZeros(h);
//         hashes[i] = index;
//         //std::cout << " hash " << h << " leading zeros " << index << std::endl;
//     }
//     return repetitions_;
// }

// void FM::Insert(long id) {
//     // compute one hash x = h(v), where we hash int to binary domain
//     for (int i = 0; i < repetitions_; i++) {
//         // compute hash
//         int h = ((hashab_[i][0] * id + hashab_[i][1]) % M_ + M_) % M_;
//         int index = trailingZeros(h);
//         r_bitmap[i][index] = true;
//     }
// }

// void FM::Insert(long id, int count) {
//     throw "FM does not support count insert";
// }

// void FM::Insert(long id, uint* hashes) {
//     FM::Insert(id);
// }

// void FM::Insert(long id, int count, uint* hashes) {
//     throw "FM does not support count insert";
// }

// int FM::QueryItem(long id) {
//     int sum = 0;
//     for (int i = 0; i < repetitions_; i++) {
//         int leftMostZero = 0;
//         for (int j = 0; j < b_; j++) {
//             if (r_bitmap[i][j] == false) {
//                 leftMostZero = j;
//                 break;
//             }
//         }
//         sum += leftMostZero;
//     }
//     return (int) pow(2, sum / repetitions_) * 1.2928;
// }

// int FM::QueryItem(long id, uint* hashes) {
//     return FM::QueryItem(id);
// }

// FM MergeFM(FM a, Sketch *b) {
//     FM merged = FM();
//     // first check if hashab is the same
//     // bool sameHashab = true;
//     // for (int i = 0; i < a.repetitions_; i++) {
//     //     if (a.hashab_[i][0] != b->hashab_[i][0] || a.hashab_[i][1] != b->hashab_[i][1]) {
//     //         sameHashab = false;
//     //         break;
//     //     }
//     // }
//     // if (!sameHashab) {
//     //     throw "Hashab not the same, FMs cannot be merged.";
//     // }
    
//     for (int i = 0; i < a.repetitions_; i++) {
//         for (int j = 0; j < a.b_; j++) {
//             merged.r_bitmap[i][j] = a.r_bitmap[i][j] || b->r_bitmap[i][j];
//         }
//     }
//     return merged;
// }
// // FM MergeFM(FM &a, Sketch &b) {
// //     FM merged = FM();
// //     for (int i = 0; i < a.repetitions_; i++) {
// //         for (int j = 0; j < a.b_; j++) {
// //             merged.r_bitmap[i][j] = a.r_bitmap[i][j] || b.r_bitmap[i][j];
// //         }
// //     }
// //     return merged;
// // }