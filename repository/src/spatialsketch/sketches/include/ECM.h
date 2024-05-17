#ifndef ECM_SKETCH_H_
#define ECM_SKETCH_H_

#include "Sketch.h"

#include <cstring>
#include <iostream>
#include <list>
#include <vector>

#define WINDOW_SIZE     300000


/*typedef struct Bucket {
	int exponent;
	int end;
	int start;
    Bucket(): exponent(-1), end(-1), start(-1) {};
} Bucket;

typedef struct ExpHist {
	Bucket bucket[100];
	int number;
    ExpHist() { memset(bucket, -1, sizeof(Bucket)*100); number = -1; };
} ExpHist;*/


typedef struct Bucket {
  int end;
  int start;
  Bucket(int end, int start): end(end), start(start) {};  
} Bucket;

typedef struct ExpHist {
  std::vector<std::list<Bucket>> buckets;
  ExpHist() {};
} ExpHist;


ExpHist MergeECM(std::vector<ExpHist> hists, int k_);
int HistSum(const ExpHist &hist, int t);

class ECM : public Sketch{
  public:
    int width_;

    ECM(float epsilon = 0.1, float delta = 0.05, long** hashab = nullptr);
    ~ECM();

    inline long** GetHashesCoeffLong() {
        return hashab_;
    }

    void Insert(long id) { throw std::runtime_error("ECM does not support Insert without timestamp"); }
    void Insert(long id, uint* hashes) { throw std::runtime_error("ECM does not support Insert without timestamp"); }
    void Insert(long id, int count, uint* hashes) { throw std::runtime_error("ECM does not support Insert without timestamp"); }
    void Insert(long id, long* hashes) { throw std::runtime_error("ECM does not support Insert without timestamp"); }
    void Insert(long id, int val, dyadic_cm_precompute* precompute) { throw std::runtime_error("ECM does not support Insert without timestamp"); }

    void Insert(long id, int count, long* hashes);
    void Insert(long id, int count) { long* hashes = nullptr; Insert(id, count, hashes); };

    int QueryItem(long id) { return QueryItem(id, 0); }
    int QueryItem(long id, uint* hashes) { return QueryItem(id, 0); }
    int Query(dyadic1D dyad) { throw std::runtime_error("ECM does not support Query(dyadic1D dyad)"); }
    int QueryItem(long id, long* hashes) { return QueryItem(id, 0, hashes); };
    int QueryItem(long id, int timestamp);
    int QueryItem(long id, int timestamp, long* hashes);

    // ALternatively, the hashes can be precomputed and passed to all other functions
    int GetItemHashes(long id, long* hashes);
    int GetItemHashes(long id, uint* hashes) { return GetItemHashes(id, (long*)hashes); }

    // Get memory
    inline float GetSize() {
        return size_;
    }

    bool isPrime(long prime);

    int** GetHashesCoeff() { throw std::runtime_error("ECM does not support GetHashesCoeff, only Longs"); }
    void PrecomputeInsert(long point, dyadic_cm_precompute *precompute) {};

//  protected:
    std::vector<ExpHist> GetHists(long id);
    std::vector<ExpHist> GetHists(long id, long* hashes);
    int GetK() { return k_; }


  private:
    long **hashab_;
    ExpHist **cm_;
    long size_ = 0;
    int k_;

    void InsertBucket(const int &i, const int &j, const int &t);
};

#endif  // ECM_SKETCH_H_