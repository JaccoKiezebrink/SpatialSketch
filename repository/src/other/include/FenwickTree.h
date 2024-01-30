#ifndef FENWICK_TREE_H_
#define FENWICK_TREE_H_

#include "CountMin.h"
#include "Statistics.h"

#include <vector>
#include <memory>


class FenwickTree {
    public:
        FenwickTree(int n);
        ~FenwickTree();

        // Updates & queries
        void Update(int x, int y, int item=0, int value=1);
        int QueryValue(int x, int y, int item=0, std::shared_ptr<statistics> stats=nullptr);
        std::vector<int> QueryCounters(int x, int y, int item=0, std::shared_ptr<statistics> stats=nullptr);
        int RangeQuery(int x1, int y1, int x2, int y2, int item=0, std::shared_ptr<statistics> stats=nullptr);

        // Get size of Fenwick tree in bytes
        float GetSize() {
            // Todo: probably incorrect
            return n_ * n_ * grid_[0][0].GetSize();
        }

    private:
        int n_;
        int d_;
        CountMin **grid_;
};

#endif  // FENWICK_TREE_H_