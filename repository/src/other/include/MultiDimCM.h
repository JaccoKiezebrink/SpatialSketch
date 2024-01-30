#ifndef MULTI_DIM_CM_H_
#define MULTI_DIM_CM_H_

#include "CountMin.h"
#include "DyadicRanges.h"
#include "Utils.h"

#include <iostream>

// Multidimensional Count-Min sketch of 3 Dimensions [x, y, ip] of ranges N x N x 2^32
class MultiDimCM {
    public:
        MultiDimCM(int N, float epsilon, float delta, long memory_limit=0);
        ~MultiDimCM();

        void Update(long lon, long lat, long ip, int value=1);
        long Query(std::vector<range> ranges, long ip_start, long ip_end);

        uint GetSize() {
            return n_dim_ * n_dim_ * ip_dim_ * cells_[0][0][0]->GetSize();
        }

        void Print();

    private:
        int n_;
        float epsilon_;
        float delta_;
        CountMin**** cells_ = nullptr;
        int n_dim_;  // log2(n_) + 1
        int ip_dim_ = 33;

        int nr_hashes_ = 0;
        uint* hashes_ = nullptr;

        inline long DimToKey(int x, int y, long ip) {
            /*long long x_ = x;
            long long y_ = y;
            return long ((x_ + y_ * (n_+1)) << 32 | ip);*/

            // (long) x + (long) y + ip;
            //(((x + y + ip) % LONG_MAX) + LONG_MAX) % LONG_MAX;
            return  (x + y * (4096)) + ip * ((long) UINT_MAX);
        }

        void ConcatVectors(std::vector<dyadic1D> &a, std::vector<dyadic1D> &b);
        std::vector<dyadic1D> ObtainIntervals(dyadic1D target, dyadic1D base);
        std::vector<dyadic3D> GetDyadicIntervals(int x1, int y1, int x2, int y2, long ip_start, long ip_end);
};

#endif  // MULTI_DIM_CM_H_