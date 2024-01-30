#include "MultiDimCM.h"

#include <cmath>
#include <iostream>

MultiDimCM::MultiDimCM(int N, float epsilon, float delta, long memory_limit) {
    n_ = N;
    epsilon_ = epsilon;
    delta_ = delta;
    n_dim_ = std::log2(n_) + 1;

    float eps_correction = 1.0f;
    if (memory_limit > 0) {
        epsilon_ = (((int) ceil(log(1 / delta)) * n_dim_ * n_dim_ * ip_dim_) / ((double) memory_limit)) * std::exp(1) * sizeof(int);
        std::cout << "Initializing MultiDimCM with epsilon " << epsilon_ << " under mem constraint " << ceil(memory_limit / 1024 / 1024) << "MB" << std::endl;
    } else {
        eps_correction = (n_dim_ - 1) * (n_dim_ - 1) * (ip_dim_ - 1);  // Every cm's epsilon has to be set to e' = e / (#levels)
    }

    cells_ = new CountMin***[n_dim_];
    for (int i = 0; i < n_dim_; i++) {
        cells_[i] = new CountMin**[n_dim_];
        for (int j = 0; j < n_dim_; j++) {
            cells_[i][j] = new CountMin*[ip_dim_];
            for (int k = 0; k < ip_dim_; k++) {
                cells_[i][j][k] = new CountMin(epsilon_ / eps_correction, delta_);
            }
        }
    }
}

MultiDimCM::~MultiDimCM() {
    for (int i = 0; i < n_dim_; i++) {
        for (int j = 0; j < n_dim_; j++) {
            delete[] cells_[i][j];
        }
        delete[] cells_[i];
    }
    delete[] cells_;
}


// Verify if target lays within interval, for example target = 5, interval = [3, 7] returns true
inline bool IsWithin(long target, std::pair<long, long> interval) {
    return target >= interval.first && target <= interval.second;
}

inline bool IsWithin(long target, long start, long end) {
    return target >= start && target <= end;
}


std::vector<std::pair<long, long>> FindChildInterval(long target, long int_start, long int_end, int levels) {
    long start = int_start;
    long end = int_end;
    std::vector<std::pair<long,long>> intervals(levels, std::make_pair(-1,-1));
    long diff = end - start + 1;
    for (int i = 0; i < levels; i++) {
        intervals[i].first = start;
        intervals[i].second = end;
        if (target == start && target == end) {
            // base case (target, target)
            break;
        } else if (start == end && target != start) {
            std::cout << "Error in base case " << target << ", (" << start << ", " << end << ")" << std::endl;
            break;
        } else {
            // Split interval
            if (IsWithin(target, start, start + diff/2 - 1)) {
                end = start + diff/2 - 1;
                //intervals[i].first = start;
                //intervals[i].second = end;
            } else if (IsWithin(target, start + diff/2, end)) {
                start = start + diff/2;
                //intervals[i].first = start;
                //intervals[i].second = end;   
            } else {
                std::cout << "Error in recursion " << target << ", (" << start << ", " << end << ")" << std::endl;
                break;
            }
            diff = diff/2;
        }
    }
    return intervals;
}

void MultiDimCM::Update(long lon, long lat, long ip, int value) {
    //nr_hashes_ = cells_[0][0][0]->GetItemHashes(DimToKey, hashes_);

    std::vector<std::pair<long, long>> x_intervals, y_intervals, ip_intervals;
    x_intervals = FindChildInterval(lon + 1, 1, n_, n_dim_);
    y_intervals = FindChildInterval(lat + 1, 1, n_, n_dim_);
    long ip_end = UINT_MAX;
    ip_intervals = FindChildInterval(ip + 1, 1, ip_end + 1, ip_dim_);

    // Combine them into 3d intervals and update them
    for (auto x_int : x_intervals) {  // go over intervals from smallest to largest
        for (auto y_int : y_intervals) { 
            for (auto ip_int : ip_intervals) { 
                //x_int.first--; x_int.second--; y_int.first--; y_int.second--; ip_int.first--; ip_int.second--;  // convert to 0-indexed
                long key = DimToKey((x_int.first-1) / (x_int.second - x_int.first + 1), (y_int.first-1) / (y_int.second - y_int.first + 1), (ip_int.first-1)  / (ip_int.second - ip_int.first + 1)); // TODO: verify
                //std::cout << x.first / (x.second - x.first + 1) << ", " << y.first / (y.second - y.first + 1) << ", " << ip.first  / (ip.second - ip.first + 1) << ", key " << key << std::endl;
                if (((int) std::log2(x_int.second - x_int.first + 1) == 0) && ((int) std::log2(y_int.second - y_int.first + 1) == 0) && ((int) std::log2(ip_int.second - ip_int.first + 1) == 0)) {
                    //std::cout << key << std::endl;
                }
                cells_[(int) std::log2(x_int.second - x_int.first + 1)][(int) std::log2(y_int.second - y_int.first + 1)][(int) std::log2(ip_int.second - ip_int.first + 1)]->Insert(key, value);
                //std::cout << "Inserted " << key << " in " << (int) std::log2(x.second - x.first + 1) << ", " << (int) std::log2(y.second - y.first + 1) << ", " << (int) std::log2(ip_int.second - ip_int.first + 1) << std::endl;
            }
        }
    }
    nr_hashes_ = 0;  // reset
}




// Concatenate vector b to vector a
void MultiDimCM::ConcatVectors(std::vector<dyadic1D> &a, std::vector<dyadic1D> &b) {
    a.reserve(a.size() + b.size());
    for (auto i : b) {
        a.push_back(i);
    }
}

// How does interval [start1, end1] overlap in [start2, end2]
inline int IntervalOverlap(long start1, long end1, long start2, long end2) {
    //std::cout << "IntervalOverlap [" << start1 << ", " << end1 << "] in [" << start2 << ", " << end2 << "]: ";
    if ((start2 >= start1 && end2 <= end1)) { 
        // Covers
        return 1;
    } else if (start1 >= start2 && end1 <= end2) {
        // Contained within
        return 2;
    } else if (start1 <= start2 && end1 >= start2 && end1 <= end2) {
        // Lower contained
        return 3;
    } else if (start1 >= start2 && start1 <= end2 && end1 >= end2) {
        // Upper contained
        return 4;
    } else {
        return 0;
    }
}

// Given a target 1D range that overlaps with a higher level dyadic interval in some manner, 
// obtain the dyadic intervals that together compose the target part that overlaps with the base
std::vector<dyadic1D> MultiDimCM::ObtainIntervals(dyadic1D target, dyadic1D base) {
    // Check for resolution compliance
    if (base.start == base.end) {
        std::vector<dyadic1D> res = {dyadic1D(base.start, base.end)};
        return res;
    } else if (target.start == base.start && target.end == base.end) {
        // If exactly overlap, return
        // split interval, if overlap lower, recurse lower, if overlap upper, recurse upper
        std::vector<dyadic1D> res = {target};
        return res;
    } else {
        // Recursion
        int power = (int) std::log2(base.end - base.start + 1);

        // Split interval
        dyadic1D lower_base = dyadic1D(base.start, base.end - (long) std::pow(2, power-1));
        dyadic1D upper_base = dyadic1D(base.start + (long) std::pow(2, power-1), base.end);

        std::vector<dyadic1D> lower_res = {}, upper_res = {};
        if (IntervalOverlap(target.start, target.end, lower_base.start, lower_base.end)) {
            lower_res = ObtainIntervals(dyadic1D(std::max(target.start, lower_base.start), std::min(target.end, lower_base.end)), lower_base);
        }
        if (IntervalOverlap(target.start, target.end, upper_base.start, upper_base.end)) {
            upper_res = ObtainIntervals(dyadic1D(std::max(target.start, upper_base.start), std::min(target.end, upper_base.end)), upper_base);
        }
        ConcatVectors(lower_res, upper_res);
        return lower_res;
    }
}


// Given a rectangular range, obtain the minimum number dyadic intervals composing it
std::vector<dyadic3D> MultiDimCM::GetDyadicIntervals(int x1, int y1, int x2, int y2, long ip_start, long ip_end) {
    std::vector<dyadic1D> x_intervals, y_intervals, ip_intervals;
    dyadic1D interval;

    // Sanity check
    if (x1 > x2 || y1 > y2 || ip_start > ip_end) {
        return std::vector<dyadic3D>();
    }

    // X 
    interval = dyadic1D(1, n_);
    x_intervals = ObtainIntervals(dyadic1D(x1+1, x2+1), interval);

    // Y
    interval = dyadic1D(1, n_);
    y_intervals = ObtainIntervals(dyadic1D(y1+1, y2+1), interval);

    // IP
    long end = UINT_MAX;
    interval = dyadic1D(1, end + 1);
    ip_intervals = ObtainIntervals(dyadic1D(ip_start+1, ip_end+1), interval);




    // Combine into 3d intervals
    std::vector<dyadic3D> dyadic_intervals;
    dyadic_intervals.reserve(x_intervals.size() * y_intervals.size() * ip_intervals.size());
    for (auto int_x : x_intervals) {
        for (auto int_y : y_intervals) {
            for (auto int_ip : ip_intervals) {
                dyadic3D r;
                r.x1 = int_x.start;
                r.y1 = int_y.start;
                r.x2 = int_x.end;
                r.y2 = int_y.end;
                r.ip_start = int_ip.start;
                r.ip_end = int_ip.end;
                //r.coverage = int_x.coverage * int_y.coverage * int_ip.coverage; // not relevant
                dyadic_intervals.push_back(r);
            }
        }
    }

    return dyadic_intervals;
}

long MultiDimCM::Query(std::vector<range> ranges, long ip_start, long ip_end) {
    long sum = 0;
    std::vector<dyadic3D> dyadic_intervals;
    dyadic_intervals.reserve(n_dim_ * n_dim_ * ip_dim_);
    nr_hashes_ = 0;

    // Query the sketch of every dyadic interval and accumulate the sum
    for (range r : ranges) {
        dyadic_intervals = GetDyadicIntervals(r.x1, r.y1, r.x2, r.y2, ip_start, ip_end);
        for (dyadic3D di : dyadic_intervals) {
            di.x1--; di.x2--; di.y1--; di.y2--; di.ip_start--; di.ip_end--;  // convert to 0-indexed
            long key = DimToKey(di.x1 / (di.x2 - di.x1 + 1), di.y1 / (di.y2 - di.y1 + 1), di.ip_start / (di.ip_end - di.ip_start + 1)); // TODO: verify
            //std::cout << di.x1 / (di.x2 - di.x1 + 1) << ", " << di.y1 / (di.y2 - di.y1 + 1) << ", " << di.ip_start / (di.ip_end - di.ip_start + 1) << ", key " << key << std::endl;
            if (di.x1 <= 1 && di.x2 >= 1 && di.y1 <= 1 && di.y2 >= 1 && di.ip_start <= 1 && di.ip_end >= 1) {
               // std::cout << "";
            }
            sum += cells_[(int) std::log2(di.x2 - di.x1 + 1)][(int) std::log2(di.y2 - di.y1 + 1)][(int) std::log2(di.ip_end - di.ip_start + 1)]->QueryItem(key);
        }
    }

    return sum;   
}

void MultiDimCM::Print() {
    for (int i = 0; i < n_dim_; i++) {
        for (int j = 0; j < n_dim_; j++) {
            for (int k = 0; k < ip_dim_; k++) {
                std::cout << "x: " << i << " y: " << j << " ip: " << k << " - ";
                cells_[i][j][k]->PrintCounters();
            }
        }
    }
}
