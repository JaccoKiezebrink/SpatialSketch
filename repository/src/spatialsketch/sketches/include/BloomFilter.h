#ifndef BLOOM_FILTER_H_
#define BLOOM_FILTER_H_

#include "Sketch.h"

#include <stdlib.h>
#include <vector>
#include <math.h>

class BloomFilter : public Sketch {
    public:
      //int domainSize; 
      long bits_;
      //int repetitions_;

      std::vector<bool> bloomFilter;

      BloomFilter(float delta=0.05, int domain_size=50000, long** hashab=nullptr);//1048576
      ~BloomFilter();

      int GetItemHashes(long id, uint* hashes);
      int GetItemHashes(long id, long* hashes) {
        return 0;
      }
      
      void Insert(long id);
      void Insert(long id, uint* hashes);
      void Insert(long id, int count);
      void Insert(long id, int count, uint* hashes);
      void Insert(long id, int count, long* hashes) { Insert(id); };

      int QueryItem(long id); // 1: true, 0: false. Not a boolean to adhere to abstract class.
      int QueryItem(long id, uint* hashes);
      int QueryItem(long id, long* hashes) { return QueryItem(id); };
      int QueryItem(long id, int timestamp) { return QueryItem(id); }
      int QueryItem(long id, int timestamp, long* hashes) { return QueryItem(id, hashes); }

      int** GetHashesCoeff() {
          throw "BF does not support GetHashesCoeff, only Longs";
          return nullptr;
      }
      long** GetHashesCoeffLong() {
        return hashab_;
      }
      int Query(dyadic1D dyad) {return 0;};

      inline float GetSize() {
          return ceil(bits_/8);
      }


        void Insert(long id, long* hashes) {
            Insert(id);
        }

      void Insert(long id, int val, dyadic_cm_precompute* precompute) {
          Insert(id, val);
      }
      void PrecomputeInsert(long point, dyadic_cm_precompute *precompute) {};
      long **hashab_;
      bool isPrime(long prime);
            
};

#endif // BLOOM_FILTER_H_