#ifndef CDARQ_H_
#define CDARQ_H_

#include "RangeTree.h"

#include <vector>
#include <random>
#include <queue>

// Determine if a point x lies within 2d range
inline bool inRange(std::vector<int> lower, std::vector<int> upper, std::vector<int> x) {
	if (x[0] >= lower[0] && x[0] <= upper[0] && x[1] >= lower[1] && x[1] <= upper[1])
		return true;
	return false;
}

typedef struct block {
    int count = 0;
    //bool core = false; 
    int nearest_core = -1; // nearest core group or the coregroup it belongs to if core = 
    
    block(int count = 0, int nearest_core = -1) : count(count), nearest_core(nearest_core) {};
} block;

typedef struct heap_entry {
    int count = 0;
    int x = 0;
    int y = 0;
    int z = 0;
    heap_entry() {};
    heap_entry(int count, int x, int y, int z) : count(count), x(x), y(y), z(z) {};
} heap_entry;

struct CompHeapEntry {
    bool operator() (heap_entry a, heap_entry b) {
        return a.count < b.count;
    }
};


// A class for Max Heap adapted from https://www.geeksforgeeks.org/binary-heap/
class MaxHeap {
	std::vector<heap_entry> harr; // pointer to array of elements in heap
	//int capacity; // maximum possible size of min heap
	int heap_size; // Current number of elements in min heap
public:
	// Constructor
	MaxHeap(int capacity);

    void PrintHeap();

    int GetSize(); // get size in bytes

    int GetHeapSize() { return heap_size; }

	// to heapify a subtree with the root at given index
	void MaxHeapify(int i);

	int parent(int i) { return (i-1)/2; }

	// to get index of left child of node at index i
	int left(int i) { return (2*i + 1); }

	// to get index of right child of node at index i
	int right(int i) { return (2*i + 2); }

	// to extract the root which is the maximum element
	heap_entry extractMax();

    // Find key in heap and return its index
    int findIndex(heap_entry k);

	// Increase key value of key at index i with one
	void increaseKey(int i);

	// Returns the minimum key (key at root) from min heap
	heap_entry getMax() { return harr[0]; }

	// Inserts a new key 'k'
	void insertKey(heap_entry k);
};

// ---------------------------------------------------------------------------------------------


typedef struct marq_item {
    float x = -1;
    float y = -1;
    float z = -1;

    marq_item(float x, float y, float z) : x(x), y(y), z(z) {}
} marq_item;

struct marq_item_comp {
  bool operator()(marq_item a, marq_item b) const { 
    if (a.x < b.x) {
        return true;
    } else if (a.x == b.x) {
        if (a.y < b.y) {
            return true;
        } else if (a.y == b.y) {
            if (a.z < b.z) {
                return true;
            }
        }
    }
    return false;
  }
};


// Reservoir sampling is an adapted version of the reservoir sampling provided by http://hadjieleftheriou.com/sketches/index.html
class ReservoirSampling {

  public:
    ReservoirSampling(unsigned long block_sample_size, int eps_1_);
    ~ReservoirSampling();

   // reservoir functions
    unsigned long getFrequency(int dim, int block, std::vector<float> lower, std::vector<float> upper);
    int getIntersectionFrequency(int x, int y, int z, std::vector<float> lower, std::vector<float> upper);
    std::vector<marq_item> getCombinedSamples(int x, int y, int z);
    int getFilteredUnionFrequency(std::vector<marq_item> unioned_samples, std::vector<float> lower, std::vector<float> upper);
    unsigned long getInputLength();
    void clear();
    void insert(int x, int y, int z, marq_item key);
    uint getSize();

  private:
    // reservoir functions
    void replace(int x, int y, int z, marq_item key);
    std::set<marq_item, marq_item_comp> getFilteredItemSet(int dim, int block, std::vector<float> lower, std::vector<float> upper);

    // Reservoir vars
    int eps_1_;
    unsigned long m_sampleSize;
    unsigned long m_t;
    unsigned long m_T;
    std::vector<std::vector<marq_item>> m_reservoir[3];
    std::uniform_real_distribution<double> unif;
    std::default_random_engine re;
};



// ---------------------------------------------------------------------------------------------

class MARQ {

  private:
    bool cdarq;  // only do cdarq
    int d_; //dimensions
    int* upper_limits_;
    int* lower_limits_;
    int eps_1_; //1\epsilon
    float theta_, theta_b_;
    int surplus_query_method_ = 0; // 0: query cold sketch once, 1: query cold sketch per cold cell
    std::shared_ptr<cdarq_info> info;

    block ***block_index_;
    int ccold_ = 0;
    int N_ = 0;
    CM_type *cold_sketch_;
    MaxHeap *cold_heap_;
    ReservoirSampling *reservoir_;

  public:
    std::vector<RangeTree::RangeTree<int,int>> rtrees;

    // Constructor
    MARQ(int dimensions, int expected_insertions, int* gupperlimits, int* glowerlimits, float epsilon=2.8f, float delta=0.5f, float theta=0.05f, int surplus_query_method=0, bool cdarq=false);
    ~MARQ();

    void PrintStatistics() {
        int node_count = 0;
        for (int i = 0; i < (int) rtrees.size(); i++) {
            node_count += info->tree_node_count[i];
        }
        std::cout << "\n--- MARQ Statistics ---\n"
                  << "#range trees: " << rtrees.size() << "\n"
                  << "with " << node_count << " nodes\n"
                  << "cold heap size: " << cold_heap_->GetHeapSize() << "\n"
                  << "hot regions: " << GetNrHotCells() << "\n"
                  << std::endl;
    }

    uint GetSize();

    int GetResolution() {
        return eps_1_;
    }

    int GetNrHotCells();

    int FindHotNeighbour(int x, int y, int z);

    // Compute block coordinates (index) from data point coordinates
    std::vector<int> fromCorToBlock(std::vector<int> cor);

    void update(std::vector<int> x);

    void updateBulk(std::vector<int> x);

    void finalizeBulk();

    void ColdToExistingTree(int core_group, heap_entry b);

    void ColdToNewTree(heap_entry b);

    void UpdateCold(int firstCor, int secondCor, int third);

    int QuerySurplus(int x1, int y1, int z1, int x2, int y2, int z2);

    long countQuery(std::vector<int>& lower, std::vector<int>& upper);

    int SampleQuery(std::vector<int>& lower, std::vector<int>& upper);

    unsigned int BlockToKey(int x, int y, int z) {
        return x * (eps_1_+1) * (eps_1_+1) + y * (eps_1_+1) + z;
    }

};


#endif  // CDARQ_H_