#include "ReservoirSampling.h"

#include <set>
#include <unordered_map>

RSampling::RSampling(float epsilon, float delta, unsigned long memory_limit) {
    if (memory_limit > 0) {
        m_sampleSize_ = floor((memory_limit - 4 * sizeof(unsigned long)) / sizeof(input_data)); 
        std::cout << "Sample size for memory limit " << memory_limit / 1024 / 1024 << "MB is " << m_sampleSize_ << std::endl << std::endl;
    } else {
        m_sampleSize_ = static_cast<unsigned long>((0.5 * std::log(2.0 / delta)) / std::pow(epsilon, 2.0));
        std::cout << "Sample size for epsilon " << epsilon << " and delta " << delta << " is " << m_sampleSize_ << std::endl << std::endl;
    }
    m_t_ = 0;
    m_T_ = 22;
    m_reservoir_ = std::vector<input_data>();
    m_reservoir_.reserve(m_sampleSize_);
    unif_ = std::uniform_real_distribution<double>(0,m_sampleSize_);  
}


RSampling::~RSampling() {

}


void RSampling::insert(input_data key) {
	if (m_t_ < m_sampleSize_) {
		m_reservoir_.push_back(key);
	} else {
		unsigned long M = std::rand() % m_sampleSize_;
		if (M < m_reservoir_.size()) {
			m_reservoir_[M] = key;
		}
	}

	m_t_++;
}

void RSampling::clear() {
    m_t_= 0;
    m_reservoir_.clear();
}

unsigned long RSampling::getInputLength() {
	return m_t_;
}

unsigned long RSampling::getFrequency(range r, long start_ip, long end_ip) {
	double f = static_cast<double>(m_t_) / static_cast<double>(m_sampleSize_);
	
	unsigned long c = 0;
	for (unsigned long i = 0; i < m_reservoir_.size(); i++) {
		if (m_reservoir_[i].ip >= start_ip && m_reservoir_[i].ip <= end_ip 
            && m_reservoir_[i].x >= r.x1 && m_reservoir_[i].x <= r.x2
            && m_reservoir_[i].y >= r.y1  && m_reservoir_[i].y <= r.y2) {
            c++;
        }
	}

	return static_cast<unsigned long>(c * f);
}

unsigned long RSampling::getCountDistinct(range r) {
	double f = static_cast<double>(m_t_) / static_cast<double>(m_sampleSize_);
	std::set<long> ips;
	for (unsigned long i = 0; i < m_reservoir_.size(); i++) {
		if (m_reservoir_[i].x >= r.x1 && m_reservoir_[i].x <= r.x2
            && m_reservoir_[i].y >= r.y1  && m_reservoir_[i].y <= r.y2) {
				ips.insert(m_reservoir_[i].ip);
			}
		}

	return static_cast<unsigned long>(f * ips.size());
}

unsigned long RSampling::getMembership(range r, long ip) {
	for (unsigned long i = 0; i < m_reservoir_.size(); i++) {
		if (m_reservoir_[i].ip == ip 
            && m_reservoir_[i].x >= r.x1 && m_reservoir_[i].x <= r.x2
            && m_reservoir_[i].y >= r.y1  && m_reservoir_[i].y <= r.y2) {
            return 1;
        }
	}

	return 0;
}

unsigned long RSampling::getL2Square(range r) {
	double f = static_cast<double>(m_t_) / static_cast<double>(m_sampleSize_);
	std::unordered_map<long, int> ips;
	std::set<long> ips_set;
	for (unsigned long i = 0; i < m_reservoir_.size(); i++) {
		if (m_reservoir_[i].x >= r.x1 && m_reservoir_[i].x <= r.x2
            && m_reservoir_[i].y >= r.y1  && m_reservoir_[i].y <= r.y2) {
			if (ips.find(m_reservoir_[i].ip) == ips.end()) {
				ips[m_reservoir_[i].ip] = 0;
			}
			ips[m_reservoir_[i].ip]++;
			ips_set.insert(m_reservoir_[i].ip);
        }
	}
	unsigned long L2_square = 0;
	for (auto it = ips.begin(); it != ips.end(); it++) {
		L2_square += it->second * it->second;
	}
	return L2_square; //static_cast<unsigned long>(f * L2_square);
}


void RSampling::replace(input_data key) {
	unsigned long M = std::rand() % m_sampleSize_;
	m_reservoir_[M] = key;
}

uint RSampling::getSize() {
	return sizeof(std::vector<input_data>) + sizeof(m_t_) + m_reservoir_.size() * sizeof(input_data);
}