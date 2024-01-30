#include "FenwickTree.h"

#include <math.h>
#include <algorithm>
#include <iostream>

FenwickTree::FenwickTree(int n) {
    // initialize grid  
    n_ = n;
    grid_ = new CountMin*[n_];
    for (int i = 0; i < n_; i++) {
        grid_[i] = new CountMin[n_];
    }
    d_ = grid_[0][0].repetitions_;
}


FenwickTree::~FenwickTree() {
    for (int i = 0; i < n_; i++) {
        delete[] grid_[i];
    }
    delete[] grid_;
}


// Query item value fenwick tree at x, y. Log statistics if pointer is given
int FenwickTree::QueryValue(int x, int y, int item, std::shared_ptr<statistics> stats) {
    int sum = 0;

    if (x < 0 || y < 0) {
        return sum;
    }

    int query_counter = 0;  // keep track of number of queries
    for(x = x + 1; x > 0; x -= x&-x) {
        for(int yy = y + 1; yy > 0; yy -= yy&-yy) {
            // Note that Fenwick tree's start from index one, thus we need to subtract 1
            sum += grid_[x - 1][yy - 1].QueryItem(item);
            query_counter++;
        }
    }

    // Log if statistics if pointer is given
    if (stats != nullptr) {
        stats->nr_of_subqueries += query_counter;
    }

    return sum;
}

// Query CM counters of item at x, y
std::vector<int> FenwickTree::QueryCounters(int x, int y, int item, std::shared_ptr<statistics> stats) {
    std::vector<int> counters(d_, 0);

    if (x < 0 || y < 0) {
        return counters;
    }
 
    int query_counter = 0;
    for (x = x + 1; x > 0; x -= x&-x) {
        for(int yy = y + 1; yy > 0; yy -= yy&-yy) {
            // Note that Fenwick tree's start from index one, thus we need to subtract 1
            std::vector<int> new_counters = grid_[x - 1][yy - 1].QueryCounters(item);
            // Add queried counters to the sum
            for (int i = 0; i < d_; i++) {
                counters[i] += new_counters[i];
            }
            query_counter++;
        }
    }

    // Log if statistics if pointer is given
    if (stats != nullptr) {
        stats->nr_of_subqueries += query_counter;
    }
    
    return counters;
}

// Range queries are done by querying the counters which are first composed and then the minimum is taken
int FenwickTree::RangeQuery(int x1, int y1, int x2, int y2, int item, std::shared_ptr<statistics> stats) {
    std::vector<int> counters = QueryCounters(x2, y2, item, stats);
    std::vector<int> add_counters, sub1_counters, sub2_counters;
    
    add_counters = QueryCounters(x1 - 1, y1 - 1, item, stats);
    sub1_counters = QueryCounters(x1 - 1, y2, item, stats);
    sub2_counters = QueryCounters(x2, y1 - 1, item, stats);

    // Combine the counters
    for (int i = 0; i < d_; i++) {
        counters[i] += add_counters[i] - sub1_counters[i] - sub2_counters[i];
    }

    return *std::min_element(counters.begin(), counters.end());
}


// Update item inserted at x, y with given value
void FenwickTree::Update(int x, int y, int item, int value) {
    for (x = x + 1; x <= n_; x += (x & -x)) {
        // Note that Fenwick tree's start from index one, thus we need to subtract 1
        for (int yy = y + 1; yy <= n_; yy += (yy & -yy)) {
            grid_[x - 1][yy - 1].Insert(item, value);
        }
    }
    return;
}