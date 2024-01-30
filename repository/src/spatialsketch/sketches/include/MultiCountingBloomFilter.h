/**************************************************/
	/*       Multi-stage Counting Bloom Filter        */
	/**************************************************/
	// This is an implementation of a multi-stage counting
	// Bloom filter, where every vector entry is a 4-bit
	// counter instead of a bit. Counters do not overflow
	// or underflow. A check is performed.
	// When SHA1 is used the domain is infinite
	// (SHA1 can hash strings of arbitrary length).
	// Nevertheless, only up to 2^16-1 counters and
	// at most 10 hash functions can be used.
	// When a Universal hash is used the domain is
	// sizeof(UniversalHash::value_type, but an arbitrary number
	// of counters and hash functions can be used.
	// The ids are converted to value_type using atoll.
	//
	// NOTICE: The usefulness of counting Bloom filters
	// with counters larger than 4 bits is "iffy".
	// For dense sets, it is possible that in order to
	// create an accurate filter with larger counters,
	// you need as much space as an exact solution would use
	// (see A. Broder and M. Mitzenmacher.
	// Network Applications of Bloom Filters: A Survey).
	// NOTICE: If you would like to use a similar sketch
	// with larger counters, please use the FastAMS sketch.
	class MultiCountingBloomFilter
	{
	public:
		MultiCountingBloomFilter(
			unsigned long counters,
			unsigned long hashes,
			HashType t = HT_UNIVERSAL
		);
		MultiCountingBloomFilter(
			unsigned long counters,
			unsigned long hashes,
			Tools::Random& r
		);
			// When using this constructor the hash type
			// of the filter is always HT_UNIVERSAL.
			// The random number generator is used to
			// produce the coefficients of the hash functions.

		MultiCountingBloomFilter(
			unsigned long counters,
			const std::vector<Tools::UniversalHash>& hashes
		);
			// When using this constructor the hash type
			// of the filter is always HT_UNIVERSAL.
			// The hash functions are explicitly provided.

		MultiCountingBloomFilter(const MultiCountingBloomFilter& in);
		MultiCountingBloomFilter(const byte* data);
		virtual ~MultiCountingBloomFilter();

		virtual MultiCountingBloomFilter& operator=(
			const MultiCountingBloomFilter& in
		);

		virtual void insert(const std::string& id, byte val = 1);
		virtual void erase(const std::string& id, byte val = 1);
		virtual void insert(
			const Tools::UniversalHash::value_type& id,
			byte val = 1
		);
		virtual void erase(
			const Tools::UniversalHash::value_type& id,
			byte val = 1
		);
		virtual void clear();
	
		virtual byte getFrequency(const std::string& id) const;
		virtual byte getFrequency(
			const Tools::UniversalHash::value_type& id
		) const;

		virtual unsigned long getVectorLength() const;
		virtual unsigned long getNumberOfHashes() const;

		virtual unsigned long getSize() const;
		virtual void getData(byte** data, unsigned long& length) const;

	private:
		HashType m_type;
		unsigned long m_counters;
		unsigned long m_hashes;
		byte* m_pFilter;
		unsigned long m_filterSize;
		std::vector<Tools::UniversalHash> m_hash;

		friend std::ostream& Sketches::operator<<(
			std::ostream& os,
			const MultiCountingBloomFilter& s
		);
		friend class CountMin;
	};

	std::ostream& operator<<(
		std::ostream& os,
		const Sketches::MultiCountingBloomFilter& s
	);

	/**************************************************/
	/*                 CountMin Sketch                */
	/**************************************************/
	// This implements the CountMin sketch as proposed in:
	// R. Motwani and P. Raghavan.
	// Randomized Algorithms
	// Cambridge International Series on Parallel Computation.
	// and
	// G. Cormode and S. Muthukrishnan
	// An Improved Data Stream Summary: The Count-Min Sketch and
	// its Applications, Journal of Algorithms 55(1), 2005
	//
	// The CountMin sketch is the same as a Multi-stage
	// Bloom filter. It provides guarantees on the returned
	// answers by enforcing specific vector sizes and number
	// of hash functions.
	// NOTICE: If you would like to use a similar sketch
	// with larger counters, please use the FastAMS sketch.
	/*class CountMin : public MultiCountingBloomFilter
	{
	public:
		virtual ~CountMin();

		CountMin(double epsilon, double delta);
		CountMin(double epsilon, double delta, Tools::Random& r);
	};
*/