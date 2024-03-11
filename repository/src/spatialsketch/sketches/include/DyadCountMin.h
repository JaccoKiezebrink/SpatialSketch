#ifndef DYAD_COUNT_MIN_H_
#define DYAD_COUNT_MIN_H_

#include "Sketch.h"
#include "CountMin.h"
#include "DyadicRanges.h"

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <string>


class DyadCountMin : public Sketch {
  public:
    int height_;
    float epsilon_;
    float delta_;

    //CountMin(int width, int depth);
    DyadCountMin(float epsilon = 0.1, float delta = 0.05, int height = 32);
    ~DyadCountMin();

    // Insert and query methods with actual hash computation
    void PrecomputeInsert(long point, dyadic_cm_precompute* precompute);
    void Insert(long id, int val);
    void Insert(long id, int val, dyadic_cm_precompute* precompute);
    
    int Query(dyadic1D dyad);

    // Get memory
    inline float GetSize() {
        std::cout << "Memory dyadiccm sketch part: " << (height_ - exact_levels_) * sketch_size_ << " - exact part: " << exact_counts_.size() * sizeof(uint) << std::endl;
        return (height_ - exact_levels_) * sketch_size_ + exact_counts_.size() * sizeof(uint);
    }

    // Functions required because this class inerhits from sketch
    void Insert(long id, int count, uint* hashes) {
      Insert(id, count);
    };

    void Insert(long id, uint* hashes) {
      Insert(id);
    };
    void Insert(long id, long* hashes) {
      Insert(id);
    };
    void Insert(long id) {
      Insert(id);
    };
    void Insert(long id, int count, long* hashes) { Insert(id, count); };
    
    int QueryItem(long id) {
      return Query((dyadic1D) id);
    };

    int QueryItem(long id, uint* hashes) {
      return Query((dyadic1D) id);
    };
    int QueryItem(long id, long* hashes) { return QueryItem(id); };
    int QueryItem(long id, int timestamp) { return QueryItem(id); }
    int QueryItem(long id, int timestamp, long* hashes) { return QueryItem(id, hashes); }

    int GetItemHashes(long id, uint* hashes) {
      return 0;
    }

    int GetItemHashes(long id, long* hashes) {
      return 0;
    }

    virtual int** GetHashesCoeff() {
      return nullptr;
    };
    virtual long** GetHashesCoeffLong(){
      return nullptr;
    }


  private:
    std::vector<CountMin*> cms_;
    std::vector<int> exact_counts_;
    int exact_levels_ = 14;  // 14 generally the minimum where exact counts size is still less than cm (for eps < 0.05)
    size_t sketch_size_ = 0;
};

#endif  // DYAD_COUNT_MIN_H_