#include "Utils.h"
#include "Postgres.h"

#include <random>
#include <iostream>


class RSampling {
    public:

        RSampling(float epsilon, float delta, unsigned long memory_limit = 0);
        ~RSampling();

        // reservoir functions
        void replace(input_data key);
        unsigned long getFrequency(range r, long start_ip, long end_ip);
        unsigned long getCountDistinct(range r);
        unsigned long getMembership(range r, long start_ip);
        unsigned long getL2Square(range r);

        unsigned long getInputLength();
        void clear();
        void insert(input_data key);
        uint getSize();
    private:
    // Reservoir vars
        unsigned long m_sampleSize_;
        unsigned long m_t_;
        unsigned long m_T_;
        std::vector<input_data> m_reservoir_;
        std::uniform_real_distribution<double> unif_;
        std::default_random_engine re_;
};


