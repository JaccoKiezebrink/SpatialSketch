#include "DyadCountMin.h"

#include <iostream>
#include <limits.h>
#include <math.h>

std::vector<std::pair<long, long>> FindChildInterval(long target);
std::vector<dyadic1D> ObtainIntervals(dyadic1D target, dyadic1D base);

DyadCountMin::DyadCountMin(float epsilon, float delta, int height) : Sketch() {
    height_ = height;
    epsilon_ = epsilon;
    delta_ = delta;

    cms_.reserve(height - exact_levels_ +1);
    int **hash_coeffs = nullptr;
    int pow = 1;
    for (int i = 0; i < height - exact_levels_ + 1; i++) {
        cms_.push_back(new CountMin(epsilon / float(height - exact_levels_), delta, hash_coeffs));
        if (hash_coeffs == nullptr) {
            hash_coeffs = cms_[0]->GetHashesCoeff();
            sketch_size_ = cms_[0]->GetSize();
        }
    }
    exact_counts_ = std::vector<int>(2*std::pow(2, exact_levels_-1)-1, 0);
}

DyadCountMin::~DyadCountMin() {
    for (int i = 0; i < (int) cms_.size() + 1; i++) {
        delete[] cms_[i];
    }
    cms_.clear();
    exact_counts_.clear();
};


void DyadCountMin::PrecomputeInsert(long point, dyadic_cm_precompute* precompute) {
    if (precompute->hashes.size() > 0) precompute->hashes.clear();
    if (precompute->intervals.size() > 0) precompute->intervals.clear();

    precompute->intervals = FindChildInterval(point);

    long pow = 1;
    long offset = 0;
    for (int i = 0; i < height_+1; i++) {
        uint* hashes = nullptr;
        if (i < exact_levels_) {
            precompute->hashes.push_back(hashes);
        } else {
            // cm key is (start interval / interval size), no offset needed as each cm is independent
            hashes = new uint[cms_[i - exact_levels_]->repetitions_];
            precompute->hashes.push_back(hashes);
            cms_[i - exact_levels_]->GetItemHashes(precompute->intervals[i].first / ((long) precompute->intervals[i].second - (long) precompute->intervals[i].first + 1), precompute->hashes[i]);
        }
    }
}


// Insert and query methods with actual hash computation
void DyadCountMin::Insert(long point, int val) {
    std::vector<std::pair<long,long>> intervals = FindChildInterval(point);
    long pow = 1;
    long offset = 0;
    for (int i = 0; i < height_+1; i++) {
        if (i < exact_levels_) {
            // Exact query computation is (start interval / interval size) + offset
            if (i == 0) {
                offset = 0;
            } else {
                offset = 2 * pow - 1;
                pow *= 2;
            }

            exact_counts_[(intervals[i].first / (intervals[i].second - intervals[i].first + 1)) + offset] += val;
        } else {
            // cm key is (start interval / interval size), no offset needed as each cm is independent
            cms_[i - exact_levels_]->Insert((intervals[i].first / ((long) intervals[i].second - (long) intervals[i].first + 1)), val);
        }
    }
}


// Insert and query methods with actual hash computation
inline void DyadCountMin::Insert(long point, int val, dyadic_cm_precompute* precompute) {
    std::vector<std::pair<long, long>> intervals = precompute->intervals;
    long pow = 1;
    long offset = 0;
    for (int i = 0; i < height_+1; i++) {
        if (i < exact_levels_) {
            // Exact query computation is (start interval / interval size) + offset
            if (i == 0) {
                offset = 0;
            } else {
                offset = 2 * pow - 1;
                pow *= 2;
            }

            exact_counts_[(intervals[i].first / (intervals[i].second - intervals[i].first + 1)) + offset] += val;
        } else {
            // cm key is (start interval / interval size), no offset needed as each cm is independent
            cms_[i - exact_levels_]->Insert((intervals[i].first / ((long) intervals[i].second - (long) intervals[i].first + 1)), val, precompute->hashes[i]);
        }
    }
}


int DyadCountMin::Query(dyadic1D dyad) {
    dyadic1D target = dyadic1D(dyad.start, dyad.end);
    dyadic1D base = dyadic1D(0, UINT32_MAX);
    //dyadic1D base = dyadic1D(std::pow(2, floor(std::log2(dyad.start))), std::pow(2, ceil(std::log2(dyad.end))) + 1);
    std::vector<dyadic1D> intervals = ObtainIntervals(target, base);
    int res = 0;
    for (int i = 0; i < (int) intervals.size(); i++) {
        int level = height_ - (int) std::log2(intervals[i].end - intervals[i].start + 1);
        if (level < exact_levels_) {
            int offset = std::pow(2, level) - 1;

            res += exact_counts_[(intervals[i].start / (intervals[i].end - intervals[i].start + 1)) + offset];
        } else {
            res += cms_[level - exact_levels_]->QueryItem((intervals[i].start / ((long) intervals[i].end - (long) intervals[i].start + 1)));
        }
    }

    return res;
}



// Verify if target lays within interval, for example target = 5, interval = [3, 7] returns true
inline bool IsWithin(long target, std::pair<long, long> interval) {
    return target >= interval.first && target <= interval.second;
}

inline bool IsWithin(long target, long start, long end) {
    return target >= start && target <= end;
}

std::vector<std::pair<long, long>> FindChildInterval(long target) {
    long start = 0;
    long end = UINT_MAX;
    std::vector<std::pair<long,long>> intervals(32 + 1, std::make_pair<long,long>(0, 0));
    long diff = (long) end - start + 1;
    for (int i = 0; i < 32 + 1; i++) {
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


// How does interval [start1, end1] overlap in [start2, end2]
int IntervalOverlap(long start1, long end1, long start2, long end2) {
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

// Concatenate vector b to vector a
void ConcatVectors(std::vector<dyadic1D> &a, std::vector<dyadic1D> &b) {
    a.reserve(a.size() + b.size());
    for (auto i : b) {
        a.push_back(i);
    }
}

std::vector<dyadic1D> ObtainIntervals(dyadic1D target, dyadic1D base) {
    // Check for resolution compliance
    if (target.start == base.start && target.end == base.end) {
        // If exactly overlap, return
        // split interval, if overlap lower, recurse lower, if overlap upper, recurse upper
        std::vector<dyadic1D> res = {target};
        return res;
    } else {
        // Recursion
        int power = (int) std::log2((long) base.end - (long) base.start + 1);

        // Split interval
        dyadic1D lower_base = dyadic1D(base.start, base.end - (long) std::pow(2, power-1));
        dyadic1D upper_base = dyadic1D(base.start + (long) std::pow(2, power-1), base.end);

        std::vector<dyadic1D> lower_res = {}, upper_res = {};
        if (IntervalOverlap(target.start, target.end, lower_base.start, lower_base.end)) {
            //std::cout << "Lower recursion with " << std::max(target.first, lower_base.first) << ", " << std::min(target.second, lower_base.second) << std::endl;
            /*if (target.end - target.start + 1 == resolution_) {
                lower_res = {target};
            } else {*/
            lower_res = ObtainIntervals(dyadic1D(std::max(target.start, lower_base.start), std::min(target.end, lower_base.end)), lower_base);
            if (lower_res.size() == 0) {
                dyadic1D partial = base;
                partial.coverage = (float) (target.end - target.start + 1) / (float) (base.end - base.start + 1);
                lower_res = {partial};
            }
        }
        if (IntervalOverlap(target.start, target.end, upper_base.start, upper_base.end)) {
            //std::cout << "Upper recursion with " << std::max(target.first, upper_base.first) << ", " << std::min(target.second, upper_base.second) << std::endl;
            //std::cout << "base " << upper_base.first << ", " << upper_base.second << std::endl;
            /*if (target.end - target.start + 1 == resolution_) {
                upper_res = {target};
            } else {*/
            upper_res = ObtainIntervals(dyadic1D(std::max(target.start, upper_base.start), std::min(target.end, upper_base.end)), upper_base);
            if (upper_res.size() == 0) {
                dyadic1D partial = base;
                partial.coverage = (float) (target.end - target.start + 1) / (float) (base.end - base.start + 1);
                upper_res = {partial}; 
            }
        }
        ConcatVectors(lower_res, upper_res);
        return lower_res;
    }
}