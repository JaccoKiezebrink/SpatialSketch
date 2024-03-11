#ifndef COUNT_MIN_H_
#define COUNT_MIN_H_

#include "Sketch.h"

#include <stdlib.h>
#include <vector>


// Add counters of b to a
inline void AddCounters(std::vector<int> &a, std::vector<int> b) {
    for (size_t i = 0; i < a.size(); i++) {
        a[i] += b[i];
    }
}

// Add counter of a to b
inline void SubtractCounters(std::vector<int> &a, std::vector<int> b) {
    for (size_t i = 0; i < a.size(); i++) {
        a[i] -= b[i];
    }
}

class CountMin : public Sketch {
  public:
    int width_;
    //int repetitions_;

    //CountMin(int width, int depth);
    CountMin(float epsilon = 0.1, float delta = 0.05, int** hashab = nullptr);
    ~CountMin();

    void PrintCounters();
    //void SetupHashes(int** hashab);
    inline int** GetHashesCoeff() {
        return hashab_;
    }

    inline long** GetHashesCoeffLong() {
        return nullptr;
    }

    // Insert and query methods with actual hash computation
    void Insert(long id);
    void Insert(long id, uint* hashes);
    void Insert(long id, int count);
    void Insert(long id, int count, uint* hashes);
    void Insert(long id, int count, long* hashes) { Insert(id, count); };

    std::vector<int> QueryCounters(long id);
    std::vector<int> QueryCounters(long id, uint* hashes);

    int QueryItem(long id);
    int QueryItem(long id, uint* hashes);
    int QueryItem(long id, long* hashes) { return QueryItem(id); };
    int QueryItem(long id, int timestamp) { return QueryItem(id); }
    int QueryItem(long id, int timestamp, long* hashes) { return QueryItem(id, hashes); }

    void Merge(Sketch *cm);
    long L2Estimate();

    // ALternatively, the hashes can be precomputed and passed to all other functions
    int GetItemHashes(long id, uint* hashes);
    int GetItemHashes(long id, long* hashes) {
        return 0;
    }

    // Get memory
    inline float GetSize() {
        return sizeof(int) * width_ * repetitions_;
    }

    int Query(dyadic1D dyad) {
        return 0;
    }

    void Insert(long id, long* hashes) {
        Insert(id);
    }

    void Insert(long id, int count, dyadic_cm_precompute* precompute) {
        Insert(id, count);
    }

    void PrecomputeInsert(long point, dyadic_cm_precompute *precompute) {
        return;
    }
    bool isPrime(long prime);


  private:
    int **hashab_;
};

#endif  // COUNT_MIN_H_