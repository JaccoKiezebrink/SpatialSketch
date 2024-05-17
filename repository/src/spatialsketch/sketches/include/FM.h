#ifndef FM_H_
#define FM_H_

#include "Sketch.h"

#include <stdlib.h>
#include <vector>
#include <math.h>


class FM : public Sketch {
    public:
        // Hyperloglog sketch
        //int repetitions_; // m: number of hash functions
        int b_=32; // bits in domain (log(M) where M is domain size)
        long M_ = 2147483647;//(2^31) - 1; // domain size
        //std::vector<std::vector<bool> > r_bitmap;


        FM(float epsilon = 0.1, float delta = 0.05, long** hashab = nullptr);
        ~FM();

        void PrintHashab();
        
        void Insert(long id);
        void Insert(long id, long* hashes);
        void Insert(long id, uint* hashes) {
            Insert(id);
        }

        void Insert(long id, int count);
        void Insert(long id, int count, uint* hashes);
        void Insert(long id, int count, long* hashes) { Insert(id, count); }

        int QueryItem(long id);
        int QueryItem(long id, uint* hashes);
        int QueryItem(long id, long* hashes) { return QueryItem(id); };
        int QueryItem(long id, int timestamp) { return QueryItem(id); }
        int QueryItem(long id, int timestamp, long* hashes) { return QueryItem(id, hashes); }

        int GetItemHashes(long id, uint* hashes) {
            return 0;
        }
        int GetItemHashes(long id, long* hashes);

        long** GetHashesCoeffLong() {
            //PrintHashab();
            return hashab_;
        }
        int** GetHashesCoeff() {
            //PrintHashab();
            throw " FM does not support GetHashesCoeff, only Longs";
            int** test = new int*[repetitions_];
            return test;
        }
        int Query(dyadic1D dyad) {return 0;}; 

        inline float GetSize() {
            return repetitions_ * ceil(b_/8);
        }

        void Insert(long id, int val, dyadic_cm_precompute* precompute) {
            Insert(id, val);
        }
        void PrecomputeInsert(long point, dyadic_cm_precompute *precompute) {};
        long** hashab_;

        bool isPrime(long prime);
        
};


void MergeFM(FM &a, Sketch *b);

#endif // FM_H_
