#include "ECM.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <limits.h>


bool ECM::isPrime(long prime) {
    if (prime <= 1) {
        return false;
    }
    if (prime <= 3) {
        return true;
    }
    if (prime % 2 == 0 || prime % 3 == 0) {
        return false;
    }
    if (prime % 5 == 0 || prime % 7 == 0) {
        return false;
    }
    if (prime % 11 == 0 || prime % 13 == 0) {
        return false;
    }
    if (prime % 17 == 0 || prime % 19 == 0) {
        return false;
    }
    return true;
}


ECM::ECM(float epsilon, float delta, long** hashab) : Sketch() {
    width_ = (int) ceil(std::exp(1) / (std::sqrt(epsilon + 1.0f) - 1.0f));
    repetitions_ = (int) ceil(log(1 / delta));
    cm_ = new ExpHist*[repetitions_];
    size_ = sizeof(int) * width_ * repetitions_; // initially all counters of the cm only consist the number of buckets
    k_ = std::ceil(1.0f/epsilon);

    if (hashab != nullptr) {
        hashab_ = hashab;
        //std::cout << "hashab passed\n";
    }

    for (int i = 0; i < repetitions_; i++) {
        cm_[i] = new ExpHist[width_];

        if (hashab == nullptr) {
            if (i == 0) {
                hashab_ = new long*[repetitions_];
            }
            hashab_[i] = new long[3];
            std::srand((i+1)*7);
            //std:: cout << std::rand() << ", " << std::rand() << ", " << std::rand() << std::endl;
            //hashab_[i][0] = int(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
            //hashab_[i][1] = int(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
            //hashab_[i][0] = int(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
            hashab_[i][0] = std::rand();//long(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
            hashab_[i][1] = std::rand();//long(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
            // generate random large prime number
            long prime = 10000 + ( std::rand() % ( RAND_MAX - 10000 + 1 ) );
            int count = 0;
            while (!isPrime(prime) && count < 1000) {
                prime = 10000 + ( std::rand() % ( RAND_MAX - 10000 + 1 ) );
            }
            
            hashab_[i][2] = prime;
        }
    }
}


ECM::~ECM() {
    for (int i = 0; i < repetitions_; i++) {
        delete[] cm_[i];
    }
    delete[] cm_;
}


// Given an item id, compute the d hashes for it to be reused
inline int ECM::GetItemHashes(long id, long* hashes) {
    for (int i = 0; i < repetitions_; i++) {
        hashes[i] = ((hashab_[i][0] * id + hashab_[i][1]) % hashab_[i][2] % width_ + width_) % width_;
    }
    return repetitions_;
}

// Count is timestamp
void ECM::Insert(long id, int count, long* hashes) {
    for (int i = 0; i < repetitions_; i++) {
        if (hashes == nullptr) {
            //ExpireBucket(i, ((hashab_[i][0] * id + hashab_[i][1]) % hashab_[i][2] % width_ + width_) % width_, timestamp);
            InsertBucket(i, ((hashab_[i][0] * id + hashab_[i][1]) % hashab_[i][2] % width_ + width_) % width_, count);
        } else {
            //ExpireBucket(i, hashes[i], timestamp);
            InsertBucket(i, hashes[i], count);
        }
    }
}



void ECM::InsertBucket(int i, int j, int t)	{
    // Add new item as singleton bucket
    Bucket new_bucket = Bucket(t, t);
    ExpHist *hist = &cm_[i][j];
    // If ExpHist is empty, intialize first exponent list
    if (hist->buckets.size() == 0) {
        std::list<Bucket> temp = std::list<Bucket>();
        temp.push_front(new_bucket);
        hist->buckets.push_back(temp);
    } else {
        // else simply add it to the front of the 1 exponent bucket
        hist->buckets[0].push_front(new_bucket);
    }
    //hist->nr_buckets++;

    // for every exponent list
    for (size_t k = 0; k < hist->buckets.size(); k++) {
        // if number of buckets exceeds k, merge the last two buckets
        if ((int) hist->buckets[k].size() > k_) {
            Bucket b1 = hist->buckets[k].back(); hist->buckets[k].pop_back();
            Bucket b2 = hist->buckets[k].back(); hist->buckets[k].pop_back();

            // If this is the last exponent list
            if ((k+1) == hist->buckets.size()) {
                hist->buckets.push_back(std::list<Bucket>());
            }
            // Create new bucket with exponent k+1
            hist->buckets[k+1].push_back(Bucket(b2.start, b1.end));
            //hist->nr_buckets--; // decrease number of buckets due to merge
        }
    }



	/*int z = cm_[i][j].number;
	int p = -1;
	int value = 0;
	int first = 0;

    // If exp hist is empty create first bucket
	if (cm_[i][j].number == -1)	{
		cm_[i][j].bucket[0].exponent = 0;
		cm_[i][j].bucket[0].start = t;
		cm_[i][j].bucket[0].end = t;
		cm_[i][j].number = 1;
        size_ += sizeof(Bucket);
	} else {
    long **hashab_;
    ExpHist **cm_;
    long size_ = 0;
		while (p < z) {
			if ((cm_[i][j].bucket[p + 2].exponent == value) && (p == -1)) {   // If neighboring bucket size can be increased, case 1
				cm_[i][j].bucket[p + 2].exponent++;
				cm_[i][j].bucket[p + 2].end = cm_[i][j].bucket[p + 1].end;
				cm_[i][j].bucket[p + 1].start = cm_[i][j].bucket[p + 2].end;
				cm_[i][j].bucket[p + 1].end = t;
				p = p + 2;
				value = cm_[i][j].bucket[p].exponent;
			} else if ((cm_[i][j].bucket[p + 2].exponent == value) && (p != -1)) {  // If neighboring bucket size can be increased, case 2
				cm_[i][j].bucket[p + 2].exponent++;
				cm_[i][j].bucket[p + 2].end = cm_[i][j].bucket[p + 1].end;
				cm_[i][j].bucket[p + 1].start = cm_[i][j].bucket[p + 2].end;
				for (int q = p + 1; q > first; q--)	{
					cm_[i][j].bucket[q].exponent = cm_[i][j].bucket[q - 1].exponent;
					cm_[i][j].bucket[q].start = cm_[i][j].bucket[q - 1].start;
					cm_[i][j].bucket[q].end = cm_[i][j].bucket[q - 1].end;
				}
				first++;
				cm_[i][j].number--;
                size_ -= sizeof(Bucket);
				p = p + 2;
				value = cm_[i][j].bucket[p].exponent;
			} else if ((cm_[i][j].bucket[p + 2].exponent != value) && (p != -1)) {  // no more merging to be done
				break;
			} else {  // shift all buckets to create a new bucket
				for (int q = z; q > 0; q--) {
					cm_[i][j].bucket[q].exponent = cm_[i][j].bucket[q - 1].exponent;
					cm_[i][j].bucket[q].start = cm_[i][j].bucket[q - 1].start;
					cm_[i][j].bucket[q].end = cm_[i][j].bucket[q - 1].end;
				}
				cm_[i][j].number++;
				cm_[i][j].bucket[0].exponent = 0;
				cm_[i][j].bucket[0].start = cm_[i][j].bucket[1].end;
				cm_[i][j].bucket[0].end = t;
                size_ += sizeof(Bucket);
				break;
			}
		}
		int qq = 0;
		if (first != 0)	{
			for (int q = first; q < first + cm_[i][j].number; ++q) {
				cm_[i][j].bucket[qq].exponent = cm_[i][j].bucket[q].exponent;
				cm_[i][j].bucket[qq].start = cm_[i][j].bucket[q].start;
				cm_[i][j].bucket[qq].end = cm_[i][j].bucket[q].end;
				qq++;
			} for (int q = qq; q < qq + first; ++q)	{
				cm_[i][j].bucket[q].exponent = -1;
				cm_[i][j].bucket[q].start = -1;
				cm_[i][j].bucket[q].end = -1;
			}
		}
	}*/
}


// Query the respective d counters for the given item id, hashes not given
std::vector<ExpHist> ECM::GetHists(long id) {
    std::vector<ExpHist> hists = std::vector<ExpHist>(repetitions_);
    for (int i = 0; i < repetitions_; i++) {
        hists[i] = cm_[i][((hashab_[i][0] * id + hashab_[i][1]) % hashab_[i][2] % width_ + width_) % width_];
    }

    return hists;
}

std::vector<ExpHist> ECM::GetHists(long id, long* hashes) {
    std::vector<ExpHist> hists = std::vector<ExpHist>(repetitions_);
    for (int i = 0; i < repetitions_; i++) {
        hists[i] = cm_[i][hashes[i]];
    }

    return hists;
}


int ECM::QueryItem(long id, int timestamp) {
    std::vector<ExpHist> hists = GetHists(id);
    int min = INT_MAX;
    for (auto ExpHist : hists) {
        int sum = HistSum(ExpHist, timestamp);
        min = std::min(min, sum);
    }
    return min;
}


int ECM::QueryItem(long id, int timestamp, long* hashes) {
    std::vector<ExpHist> hists = GetHists(id, hashes);
    int min = INT_MAX;
    for (auto ExpHist : hists) {
        int sum = HistSum(ExpHist, timestamp);
        min = std::min(min, sum);
    }
    return min;
}


int HistSum(const ExpHist &hist, int t) {
    if (hist.buckets.size() == 0) {
        return 0;
    }

    int sum = 0;
    for (int exp = 0; exp < (int) hist.buckets.size(); exp++) {  // loop over exponent lists
        // if the last bucket falls within the window, we can simply aggregate the entire list
        if (hist.buckets[exp].back().start >= t) { 
            sum += std::pow(2, exp) * hist.buckets[exp].size();
            continue;
        }
        // If we already exceeded the time window, we can break
        if (hist.buckets[exp].front().end < t) {
            break;
        }

        // otherwise, loop over buckets in list
        for (Bucket bucket : hist.buckets[exp]) {
            if (bucket.start >= t) {
                sum += std::pow(2, exp);
            } else if (bucket.start < t && bucket.end >= t) {
                sum += (std::pow(2, exp) / 2);
                break;
            }
        }
    }
    return sum;
}

void InsertExpHist(ExpHist &hist, int k_, int t) {
    // Add new item as singleton bucket
    Bucket new_bucket = Bucket(t, t);
    // If ExpHist is empty, intialize first exponent list
    if (hist.buckets.size() == 0) {
        std::list<Bucket> temp = std::list<Bucket>();
        temp.push_front(new_bucket);
        hist.buckets.push_back(temp);
    } else {
        // else simply add it to the front of the 1 exponent bucket
        hist.buckets[0].push_front(new_bucket);
    }
    //hist->nr_buckets++;

    // for every exponent list
    for (size_t k = 0; k < hist.buckets.size(); k++) {
        // if number of buckets exceeds k, merge the last two buckets
        if ((int) hist.buckets[k].size() > k_) {
            Bucket b1 = hist.buckets[k].back(); hist.buckets[k].pop_back();
            Bucket b2 = hist.buckets[k].back(); hist.buckets[k].pop_back();

            // If this is the last exponent list
            if ((k+1) == hist.buckets.size()) {
                hist.buckets.push_back(std::list<Bucket>());
            }
            // Create new bucket with exponent k+1
            hist.buckets[k+1].push_back(Bucket(b2.start, b1.end));
            //hist->nr_buckets--; // decrease number of buckets due to merge
        }
    }
}

ExpHist MergeECM(std::vector<ExpHist> hists, int k_) {
    ExpHist merged = ExpHist();

    std::vector<std::pair<int, int>> arrivals;  // arrivals consisting of <time, amount>
    for (size_t i = 0; i < hists.size(); i++) {
        for (size_t j = 0; j < hists[i].buckets.size(); j++) {
            for (Bucket b : hists[i].buckets[j]) {
                // Assume half of bucket arrived at start and half at end
                if (j != 0) {
                    arrivals.push_back(std::make_pair(b.start, std::pow(2, j - 1)));
                    arrivals.push_back(std::make_pair(b.end, std::pow(2, j - 1)));
                } else {
                    // If bucket size is one case
                    arrivals.push_back(std::make_pair(b.start, 1));
                }
            }
        }
    }

    // Sort arrivals based on time
    std::sort(arrivals.begin(), arrivals.end(), [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
        return a.first < b.first;
    });

    // Insert arrivals into one exp hist
    for (std::pair<int, int> arrival : arrivals) {
        for (int i = 0; i < arrival.second; i++) {
            InsertExpHist(merged, k_, arrival.first);
        }
    }

    return merged;
}