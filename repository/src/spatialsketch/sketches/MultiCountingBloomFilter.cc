// Sketches Library
//
// Copyright (C) 2005 Marios Hadjieleftheriou
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  Email:
//    mhadji@gmail.com


MultiCountingBloomFilter::MultiCountingBloomFilter(
	unsigned long counters, unsigned long hashes, HashType t
) : m_type(t), m_counters(counters), m_hashes(hashes)
{
	assert(CHAR_BIT % 4 == 0);

	if (m_hashes == 0)
		throw IllegalArgumentException("MultiCountingBloomFilter: number of hashes must be larger than zero.");

	if (m_counters == 0)
		throw IllegalArgumentException("MultiCountingBloomFilter: vector size must be larger than zero.");

	if (m_type == HT_UNIVERSAL)
	{
		Tools::Random r;
		for (unsigned long i = 0; i < m_hashes; i++)
			m_hash.push_back(Tools::UniversalHash(r));
	}
	else if (m_type == HT_SHA1)
	{
		if (m_counters >= std::pow(2.0, 16.0))
			throw Tools::IllegalArgumentException("MultiCountingBloomFilter: to use the SHA1 hash, vector size must be in [1, 2^16 - 1].");

		if (m_hashes > 10)
			throw Tools::IllegalArgumentException("MultiCountingBloomFilter: to use the SHA1 hash, number of hashes must be in [1, 10].");
	}
	else
	{
		throw Tools::NotSupportedException(
			"MultiCountingBloomFilter: This hash type is not supported yet."
		);
	}

	m_filterSize = static_cast<unsigned long>(std::ceil(m_counters * 4.0 / static_cast<double>(CHAR_BIT)));
	m_pFilter = new byte[m_filterSize * m_hashes];
	bzero(m_pFilter, m_filterSize * m_hashes * sizeof(byte));
}

MultiCountingBloomFilter::MultiCountingBloomFilter(
	unsigned long counters, unsigned long hashes, Tools::Random& r
) : m_type(HT_UNIVERSAL), m_counters(counters), m_hashes(hashes)
{
	assert(CHAR_BIT % 4 == 0);

	if (m_hashes == 0)
		throw Tools::IllegalArgumentException("MultiCountingBloomFilter: number of hashes must be larger than zero.");

	if (m_counters == 0)
		throw Tools::IllegalArgumentException("MultiCountingBloomFilter: vector size must be larger than zero.");

	m_filterSize = static_cast<unsigned long>(std::ceil(m_counters * 4.0 / static_cast<double>(CHAR_BIT)));
	m_pFilter = new byte[m_filterSize * m_hashes];
	bzero(m_pFilter, m_filterSize * m_hashes * sizeof(byte));

	for (unsigned long i = 0; i < m_hashes; i++)
		m_hash.push_back(Tools::UniversalHash(r));
}

MultiCountingBloomFilter::MultiCountingBloomFilter(
	unsigned long counters,
	const std::vector<Tools::UniversalHash>& hashes
) : m_type(HT_UNIVERSAL), m_counters(counters), m_hashes(hashes.size()),
    m_hash(hashes)
{
	assert(CHAR_BIT % 4 == 0);

	if (m_hashes == 0)
		throw Tools::IllegalArgumentException("MultiCountingBloomFilter: number of hashes must be larger than zero.");

	if (m_counters == 0)
		throw Tools::IllegalArgumentException("MultiCountingBloomFilter: vector size must be larger than zero.");

	m_filterSize = static_cast<unsigned long>(std::ceil(m_counters * 4.0 / static_cast<double>(CHAR_BIT)));
	m_pFilter = new byte[m_filterSize * m_hashes];
	bzero(m_pFilter, m_filterSize * m_hashes * sizeof(byte));
}

MultiCountingBloomFilter::MultiCountingBloomFilter(
	const Sketches::MultiCountingBloomFilter& in
) : m_type(in.m_type), m_counters(in.m_counters), m_hashes(in.m_hashes),
    m_filterSize(in.m_filterSize), m_hash(in.m_hash)
{
	assert(CHAR_BIT % 4 == 0);

	m_pFilter = new byte[m_filterSize * m_hashes];
	memcpy(m_pFilter, in.m_pFilter, m_filterSize * m_hashes);
}

MultiCountingBloomFilter::MultiCountingBloomFilter(const byte* data)
{
	memcpy(&m_type, data, sizeof(HashType));
	data += sizeof(HashType);
	memcpy(&m_hashes, data, sizeof(unsigned long));
	data += sizeof(unsigned long);
	memcpy(&m_counters, data, sizeof(unsigned long));
	data += sizeof(unsigned long);

	if (m_type == HT_UNIVERSAL)
	{
		for (unsigned long i = 0; i < m_hashes; i++)
		{
			Tools::UniversalHash h(data);
			m_hash.push_back(h);
			data += h.getSize();
		}
	}

	m_filterSize = static_cast<unsigned long>(std::ceil(m_counters * 4.0 / static_cast<double>(CHAR_BIT)));
	m_pFilter = new byte[m_filterSize * m_hashes];
	memcpy(m_pFilter, data, m_filterSize * m_hashes);
}

MultiCountingBloomFilter&
MultiCountingBloomFilter::operator=(
	const MultiCountingBloomFilter& in
)
{
	if (this != &in)
	{
		m_type = in.m_type;
		m_hash = in.m_hash;

		if (m_hashes != in.m_hashes || m_filterSize != in.m_filterSize)
		{
			m_counters = in.m_counters;
			m_filterSize = in.m_filterSize;
			m_hashes = in.m_hashes;
			delete[] m_pFilter;
			m_pFilter = new byte[m_filterSize * m_hashes];
		}

		memcpy(m_pFilter, in.m_pFilter, m_filterSize * m_hashes);
	}

	return *this;
}

MultiCountingBloomFilter::~MultiCountingBloomFilter()
{
	delete[] m_pFilter;
}

void MultiCountingBloomFilter::insert(const std::string& id, byte val)
{
	if (m_type == HT_UNIVERSAL)
	{
		unsigned long long l = atoll(id.c_str());
		insert(l);
	}
	else if (m_type == HT_SHA1)
	{
		Tools::SHA1Hash sha;
		unsigned long len;
		byte* data;
		sha.hash(id, &data, len);

		for (unsigned long i = 0; i < m_hashes; i++)
		{
			unsigned long h =
				*(reinterpret_cast<uint16_t*>(data + 2 * i)) % m_counters;
			unsigned long b = i * m_filterSize + (h / (CHAR_BIT / 4));
			unsigned long offset = h % (CHAR_BIT / 4);
			byte mask = (0x0F << (offset * 4));
			byte v = ((m_pFilter[b] & mask) + (val << (offset * 4))) & mask;
			m_pFilter[b] &= ~mask;
			m_pFilter[b] |= v;
		}
		delete[] data;
	}
	else
	{
		throw Tools::NotSupportedException(
			"MultiCountingBloomFilter: This hash type is not supported yet."
		);
	}
}

void MultiCountingBloomFilter::insert(
	const Tools::UniversalHash::value_type& id,
	byte val
)
{
	if (m_type == HT_UNIVERSAL)
	{
		for (unsigned long i = 0; i < m_hashes; i++)
		{
			unsigned long h = m_hash[i].hash(id) % m_counters;
			unsigned long b = i * m_filterSize + (h / (CHAR_BIT / 4));
			unsigned long offset = h % (CHAR_BIT / 4);
			byte mask = (0x0F << (offset * 4));
			byte v = ((m_pFilter[b] & mask) + (val << (offset * 4))) & mask;
			m_pFilter[b] &= ~mask;
			m_pFilter[b] |= v;
		}
	}
	else if (m_type == HT_SHA1)
	{
		std::ostringstream ss;
		ss << id << std::flush;
		insert(ss.str());
	}
	else
	{
		throw Tools::NotSupportedException(
			"MultiCountingBloomFilter: This hash type is not supported yet."
		);
	}
}

void Sketches::MultiCountingBloomFilter::erase(const std::string& id, byte val)
{
	if (m_type == HT_UNIVERSAL)
	{
		unsigned long long l = atoll(id.c_str());
		erase(l);
	}
	else if (m_type == HT_SHA1)
	{
		Tools::SHA1Hash sha;
		unsigned long len;
		byte* data;
		sha.hash(id, &data, len);

		for (unsigned long i = 0; i < m_hashes; i++)
		{
			unsigned long h =
				*(reinterpret_cast<uint16_t*>(data + 2 * i)) % m_counters;
			unsigned long b = i * m_filterSize + (h / (CHAR_BIT / 4));
			unsigned long offset = h % (CHAR_BIT / 4);
			byte mask = (0x0F << (offset * 4));
			byte v = ((m_pFilter[b] & mask) - (val << (offset * 4))) & mask;
			m_pFilter[b] &= ~mask;
			m_pFilter[b] |= v;
		}
		delete[] data;
	}
	else
	{
		throw Tools::NotSupportedException(
			"MultiCountingBloomFilter: This hash type is not supported yet."
		);
	}
}

void Sketches::MultiCountingBloomFilter::erase(
	const Tools::UniversalHash::value_type& id,
	byte val
)
{
	if (m_type == HT_UNIVERSAL)
	{
		for (unsigned long i = 0; i < m_hashes; i++)
		{
			unsigned long h = m_hash[i].hash(id) % m_counters;
			unsigned long b = i * m_filterSize + (h / (CHAR_BIT / 4));
			unsigned long offset = h % (CHAR_BIT / 4);
			byte mask = (0x0F << (offset * 4));
			byte v = ((m_pFilter[b] & mask) - (val << (offset * 4))) & mask;
			m_pFilter[b] &= ~mask;
			m_pFilter[b] |= v;
		}
	}
	else if (m_type == HT_SHA1)
	{
		std::ostringstream ss;
		ss << id << std::flush;
		erase(ss.str());
	}
	else
	{
		throw Tools::NotSupportedException(
			"MultiCountingBloomFilter: This hash type is not supported yet."
		);
	}
}

void Sketches::MultiCountingBloomFilter::clear()
{
	bzero(m_pFilter, m_filterSize * m_hashes * sizeof(byte));
}

byte Sketches::MultiCountingBloomFilter::getFrequency(
	const std::string& id
) const
{
	if (m_type == HT_UNIVERSAL)
	{
		unsigned long long l = atoll(id.c_str());
		return getFrequency(l);
	}
	else if (m_type == HT_SHA1)
	{
		byte min = std::numeric_limits<byte>::max();

		Tools::SHA1Hash sha;
		unsigned long len;
		byte* data;
		sha.hash(id, &data, len);

		for (unsigned long i = 0; i < m_hashes; i++)
		{
			unsigned long h =
				*(reinterpret_cast<uint16_t*>(data + 2 * i)) % m_counters;
			unsigned long b = i * m_filterSize + (h / (CHAR_BIT / 4));
			unsigned long offset = h % (CHAR_BIT / 4);
			byte mask = (0x0F << (offset * 4));
			byte v = ((m_pFilter[b] & mask) >> (offset * 4));
			if (v < min) min = v;
		}
		delete[] data;
		return min;
	}
	else
	{
		throw Tools::NotSupportedException(
			"MultiCountingBloomFilter: This hash type is not supported yet."
		);
	}
}

byte Sketches::MultiCountingBloomFilter::getFrequency(
	const Tools::UniversalHash::value_type& id
) const
{
	if (m_type == HT_UNIVERSAL)
	{
		byte min = std::numeric_limits<byte>::max();

		for (unsigned long i = 0; i < m_hashes; i++)
		{
			unsigned long h = m_hash[i].hash(id) % m_counters;
			unsigned long b = i * m_filterSize + (h / (CHAR_BIT / 4));
			unsigned long offset = h % (CHAR_BIT / 4);
			byte mask = (0x0F << (offset * 4));
			byte v = ((m_pFilter[b] & mask) >> (offset * 4));
			if (v < min) min = v;
		}
		return min;
	}
	else if (m_type == HT_SHA1)
	{
		std::ostringstream ss;
		ss << id << std::flush;
		return getFrequency(ss.str());
	}
	else
	{
		throw Tools::NotSupportedException(
			"MultiCountingBloomFilter: This hash type is not supported yet."
		);
	}
}

unsigned long Sketches::MultiCountingBloomFilter::getVectorLength() const
{
	return m_counters;
}

unsigned long Sketches::MultiCountingBloomFilter::getNumberOfHashes() const
{
	return m_hashes;
}

unsigned long Sketches::MultiCountingBloomFilter::getSize() const
{
	unsigned long ret =
		sizeof(HashType) +
		2 * sizeof(unsigned long) +
		m_filterSize * m_hashes * sizeof(byte);

	if (m_type == HT_UNIVERSAL)
		for (unsigned long i = 0; i < m_hash.size(); i++)
			ret += m_hash[i].getSize();

	return ret;
}

void Sketches::MultiCountingBloomFilter::getData(
	byte** data, unsigned long& length
) const
{
	length = getSize();
	*data = new byte[length];
	byte* p = *data;

	memcpy(p, &m_type, sizeof(HashType));
	p += sizeof(HashType);
	memcpy(p, &m_hashes, sizeof(unsigned long));
	p += sizeof(unsigned long);
	memcpy(p, &m_counters, sizeof(unsigned long));
	p += sizeof(unsigned long);

	unsigned long l;
	byte* buf;
	for (unsigned long i = 0; i < m_hash.size(); i++)
	{
		m_hash[i].getData(&buf, l);
		memcpy(p, buf, l);
		p += l;
		delete[] buf;
	}

	memcpy(p, m_pFilter, m_filterSize * m_hashes);
	p += m_filterSize * m_hashes;

	assert(p == (*data) + length);
}

std::ostream& Sketches::operator<<(
	std::ostream& os,
	const Sketches::MultiCountingBloomFilter& s
)
{
	os << s.m_type << " " << s.m_hashes << " " << s.m_counters;

	for (unsigned long i = 0; i < s.m_hash.size(); i++)
		os << " " << s.m_hash[i];

	for (unsigned long i = 0; i < s.m_filterSize; i++)
	{
		for (unsigned long offset = 0; offset < (CHAR_BIT / 4); offset++)
		{
			byte mask = (0x0F << (offset * 4));
			byte v = ((s.m_pFilter[i] & mask) >> (offset * 4));
			os << " " << static_cast<unsigned long>(v);
		}
	}

	return os;
}

Sketches::CountMin::~CountMin()
{
}

Sketches::CountMin::CountMin(double epsilon, double delta)
 : MultiCountingBloomFilter(static_cast<unsigned long>(std::ceil(M_E / epsilon)), static_cast<unsigned long>(std::ceil(std::log(1.0 / delta))))
{
}

Sketches::CountMin::CountMin(double epsilon, double delta, Tools::Random& r)
 : MultiCountingBloomFilter(static_cast<unsigned long>(std::ceil(M_E / epsilon)), static_cast<unsigned long>(std::ceil(std::log(1.0 / delta))), r)
{
}

