#include "marq.h"

#include "Utils.h"


// int a, b to vector
auto f = [](int a, int b) { std::vector<int> c = {a, b}; return c;};
auto f3 = [](int a, int b, int c) { std::vector<int> d = {a, b, c}; return d;};
auto f1 = [](int a) { std::vector<int> b = {a}; return b;};

// bool a, b to vector
auto g = [](bool a, bool b) { std::vector<bool> c = {a, b}; return c;};


// Prototype of a utility function to swap two integers
// A utility function to swap two elements
void swap(heap_entry *x, heap_entry *y) {
	heap_entry temp = *x;
	*x = *y;
	*y = temp;
}

// A class for Max Heap adapted from https://www.geeksforgeeks.org/binary-heap/
// Constructor: Builds a heap from a given array a[] of given size
MaxHeap::MaxHeap(int cap) {
	heap_size = 0;
	//capacity = cap;
	if (cap > (int) harr.max_size()) {
		harr.reserve(cap);
	} else {
		harr.reserve(cap);
	}
}

int MaxHeap::GetSize() {
	return sizeof(heap_entry) * heap_size + sizeof(int) * 2;
}

void MaxHeap::PrintHeap() {
	for (int i = 0; i < heap_size; i++) {
		std::cout << "[" << harr[i].count << ", " << harr[i].x << ", " << harr[i].y << ", " << harr[i].z << "]" << std::endl;
	}
	std::cout << std::endl;
}

// Inserts a new key 'k'
void MaxHeap::insertKey(heap_entry k) {
	/*if (heap_size == capacity) 	{
		std::cout << "\nOverflow: Could not insertKey [" << k.x << ", " << k.y << ", " << k.z << "]\n";
		PrintHeap();
		return;
	}*/

	// First insert the new key at the end
	heap_size++;
	int i = heap_size - 1;
	if (heap_size > (int) harr.size()) {
		harr.push_back(k);
	} else {
		harr[i] = k;
	}

	// Fix the min heap property if it is violated
	while (i != 0 && harr[parent(i)].count < harr[i].count) {
	swap(&harr[i], &harr[parent(i)]);
	i = parent(i);
	}
}

int MaxHeap::findIndex(heap_entry k) {
	for (int i = 0; i < heap_size; i++) {
		if (harr[i].x == k.x && harr[i].y == k.y && harr[i].z == k.z) {
			return i;
		}
	}
	return -1;
}

// Increase value of key at index 'i' to new_val. It is assumed that
// new_val is larger than harr[i].
void MaxHeap::increaseKey(int i) {
	harr[i].count++;
	while (i != 0 && harr[parent(i)].count < harr[i].count)	{
		swap(&harr[i], &harr[parent(i)]);
		i = parent(i);
	}
}

// Method to remove minimum element (or root) from min heap
heap_entry MaxHeap::extractMax() {
	if (heap_size <= 0)
		return heap_entry(-1, -1, -1, -1);
	if (heap_size == 1)	{
		heap_size--;
		return harr[0];
	}

	// Store the minimum value, and remove it from heap
	heap_entry root = harr[0];
	harr[0] = harr[heap_size-1];
	heap_size--;
	MaxHeapify(0);

	return root;
}


// A recursive method to heapify a subtree with the root at given index
// This method assumes that the subtrees are already heapified
void MaxHeap::MaxHeapify(int i) {
	int l = left(i);
	int r = right(i);
	int smallest = i;
	if (l < heap_size && harr[l].count > harr[i].count)
		smallest = l;
	if (r < heap_size && harr[r].count > harr[smallest].count)
		smallest = r;
	if (smallest != i) {
		swap(&harr[i], &harr[smallest]);
		MaxHeapify(smallest);
	}
}


// ------------------------------------------------------------
// Reservoir sampling is an adapted version of the reservoir sampling provided by http://hadjieleftheriou.com/sketches/index.html

ReservoirSampling::ReservoirSampling(unsigned long block_sample_size, int eps_1) {
	m_sampleSize = block_sample_size;
	eps_1_ = eps_1;
	for (int i = 0; i < 3; i++) {
		m_reservoir[i] = std::vector<std::vector<marq_item>>(eps_1_);
		for (int j = 0; j < eps_1_; j++) {
			m_reservoir[i][j].reserve(m_sampleSize);
		}
	}
	
	m_t = 0;
    m_T = 22;
    unif = std::uniform_real_distribution<double>(0,m_sampleSize);    
}

ReservoirSampling::~ReservoirSampling() {
	clear();
}

void ReservoirSampling::insert(int x, int y, int z, marq_item key) {
    static double W = std::exp(	-std::log(unif(re)) / static_cast<double>(m_sampleSize));

	static bool bInsert = false;
	static unsigned long G = 0;

	if (bInsert) {
		replace(x, y, z, key);
		bInsert = false;
	} else if (m_t < m_sampleSize) {
		m_reservoir[0][x].push_back(key);
		m_reservoir[1][y].push_back(key);
		m_reservoir[2][z].push_back(key);
	} else if (G > 0) {
		G--;
		if (G == 0) bInsert = true;
	} else if (m_t <= m_T * m_sampleSize) {
		// Algorithm X
		double V = unif(re);
		double t = m_t + 1;
		double quot = static_cast<double>(t) / static_cast<double>(t - m_sampleSize);

		while (quot > V) {
			G++;
			t++;
			quot *=
				static_cast<double>(t - m_sampleSize) / static_cast<double>(t);
		}

		if (G == 0) replace(x, y, z, key);
		else G--;
	} else {
		// Algorithm Z
		while (true) {
			unsigned long term = m_t - m_sampleSize + 1;
			double U = (double) unif(re);
			double X = static_cast<double>(m_t) * (W - 1.0);
			G = static_cast<unsigned long>(std::floor(X));
			double lhs =
				std::exp(std::log(((U * std::pow((m_t + 1.0) /
				static_cast<double>(term), 2.0)) * (term + G)) /
				(m_t + X)) / static_cast<double>(m_sampleSize));
			double rhs =
				(((m_t + X) / (static_cast<double>(term + G))) * term) /
				static_cast<double>(m_t);

			if (lhs < rhs) {
				W = rhs / lhs;
				break;
			}
		
			double y = (((U * (m_t + 1)) /	static_cast<double>(term)) * (m_t + G + 1.0)) / (m_t + X);
			unsigned long denom, numer_lim;

			if (m_sampleSize < G) {
				denom = m_t;
				numer_lim = term + G;
			} else {
				denom = m_t - m_sampleSize + G;
				numer_lim = m_t + 1;
			}
			
			for (unsigned long numer = m_t + G; numer >= numer_lim; numer--) {
				y = (y * numer) / static_cast<double>(denom);
				denom--;
			}
	
			W =	std::exp(-std::log(unif(re)) / static_cast<double>(m_sampleSize));

			if (std::exp(std::log(y) / static_cast<double>(m_sampleSize)) <= (m_t + X) / static_cast<double>(m_t)) break;
		}

		if (G == 0) replace(x, y, z, key);
		else G--;
	}

	m_t++;
}

void ReservoirSampling::clear() {
    m_t = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < eps_1_; j++) {
			m_reservoir[i][j].clear();
		}
	}
}

unsigned long ReservoirSampling::getInputLength() {
	return m_t;
}

unsigned long ReservoirSampling::getFrequency(int dim, int block, std::vector<float> lower, std::vector<float> upper) {
	double f = static_cast<double>(m_t) / static_cast<double>(m_sampleSize);
	
	unsigned long c = 0;
	for (auto sample : m_reservoir[dim][block]) {
		if (sample.x >= lower[0] && sample.x <= upper[0]
            && sample.y >= lower[1] && sample.y <= upper[1]
			&& sample.z >= lower[2] && sample.z <= upper[2] ) {
            c++;
        }
	}

	return static_cast<unsigned long>(c * f);
}


std::set<marq_item, marq_item_comp> ReservoirSampling::getFilteredItemSet(int dim, int block, std::vector<float> lower, std::vector<float> upper) {
	std::set<marq_item, marq_item_comp> result;
	for (auto sample : m_reservoir[dim][block]) {
		if (sample.x >= lower[0] && sample.x <= upper[0]
			&& sample.y >= lower[1] && sample.y <= upper[1]
			&& sample.z >= lower[2] && sample.z <= upper[2] ) {
			result.insert(sample);
		}
	}
	return result;
}

std::vector<marq_item> ReservoirSampling::getCombinedSamples(int x, int y, int z) {
	std::vector<marq_item> combined_samples = m_reservoir[0][x];
	combined_samples.insert(combined_samples.end(), m_reservoir[1][y].begin(), m_reservoir[1][y].end());
	combined_samples.insert(combined_samples.end(), m_reservoir[2][z].begin(), m_reservoir[2][z].end());

	return combined_samples;
}

int ReservoirSampling::getFilteredUnionFrequency(std::vector<marq_item> unioned_samples, std::vector<float> lower, std::vector<float> upper) {
	int count = 0;
	for (auto sample : unioned_samples) {
		if (sample.x >= lower[0] && sample.x <= upper[0]
			&& sample.y >= lower[1] && sample.y <= upper[1]
			&& sample.z >= lower[2] && sample.z <= upper[2] ) {
			count++;
		}
	}
	return count;
}

int ReservoirSampling::getIntersectionFrequency(int x, int y, int z, std::vector<float> lower, std::vector<float> upper) {
	double f = static_cast<double>(m_t) / static_cast<double>(m_sampleSize);
	
	unsigned long c = 0;
	std::set<marq_item, marq_item_comp> x_set = getFilteredItemSet(0, x, lower, upper);
	std::set<marq_item, marq_item_comp> y_set = getFilteredItemSet(1, y, lower, upper);
	std::set<marq_item, marq_item_comp> z_set = getFilteredItemSet(2, z, lower, upper);

	std::set<marq_item, marq_item_comp> xy_intersect, xyz_intersect;
	std::set_intersection(x_set.begin(), x_set.end(), y_set.begin(), y_set.end(), std::inserter(xy_intersect, xy_intersect.begin()), marq_item_comp());
	std::set_intersection(xy_intersect.begin(), xy_intersect.end(), z_set.begin(), z_set.end(), std::inserter(xyz_intersect, xyz_intersect.begin()), marq_item_comp());
	
	return (int) xyz_intersect.size();
}

void ReservoirSampling::replace(int x, int y, int z, marq_item key) {
	unsigned long M = std::rand() % m_sampleSize;
	m_reservoir[0][x][M] = key;
	m_reservoir[1][y][M] = key;
	m_reservoir[2][z][M] = key;
}

uint ReservoirSampling::getSize() {
	size_t size = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < eps_1_; j++) {
			size += m_reservoir[i][j].size() * sizeof(marq_item);
		}
	}
	return 4 * sizeof(unsigned long) + size;
}


// ------------------------------------------------------------



MARQ::MARQ(int dimensions, int expected_insertions, int* gupperlimits, int* glowerlimits, float epsilon, float delta, float theta, int surplus_query_method, bool cdarq) {
	// Initialize the sketch based on user-supplied size
	d_ = dimensions;
	eps_1_ = (int) std::ceil(1/(epsilon/2));  // marq uses eps/2 instance of cdarq
	upper_limits_ = gupperlimits;
	lower_limits_ = glowerlimits;
	theta_ = theta;	
	theta_b_ = std::min(theta, theta * (float) std::pow(epsilon / 2, d_) * (float) expected_insertions);
	surplus_query_method_ = surplus_query_method;

	if (dimensions != 3) {
		std::cout << "Only 3D is supported" << std::endl;
		exit(1);
	}

	block_index_ = new block**[eps_1_];
	for (int i = 0; i < eps_1_; i++) {
		block_index_[i] = new block*[eps_1_];
		for (int j = 0; j < eps_1_; j++) {
			block_index_[i][j] = new block[eps_1_];
		}
	}

	// Initialize cold heap and sketch
	cold_heap_ = new MaxHeap((int) std::ceil((1-theta_) * std::pow(eps_1_, d_))); // todo, (int) std::ceil((1-theta_b_) * std::pow(eps_1_, d_))
	cold_sketch_ = CM_Init((int) std::ceil(std::exp(1) / (epsilon/2)), (int) std::ceil(std::log(1 / delta)), 32);

	// Tree sketch and relevant info passed to all trees and their nodes
	info = std::make_shared<cdarq_info>();
	float e1 = (epsilon/2) / std::log(theta_b_ * (1/(epsilon/2)));  // set e1 according to paper
	info->sketch = CM_Init(std::max(1, (int) (std::ceil(std::exp(1) / e1))), std::max(1, (int) std::ceil(std::log(1 / delta))), 32);
	info->tree_count = 0;

	// given epsilon, delta, ss is x, set to eps/4d and resulting sample size divided over eps_1_ per index
	unsigned long block_sample_size = static_cast<unsigned long>(std::ceil(((0.5 * std::log(2.0 / delta)) / std::pow(epsilon / 4*dimensions, 2.0)) / (float) eps_1_));
	reservoir_ = new ReservoirSampling(block_sample_size, eps_1_);
}


MARQ::~MARQ() {
	for (int i = 0; i < eps_1_; i++) {
		for (int j = 0; j < eps_1_; j++) {
			delete[] block_index_[i][j];
		}
		delete[] block_index_[i];
	}
	delete[] block_index_;
	delete cold_heap_;
	delete reservoir_;
}


uint MARQ::GetSize() {
	uint tree_sketch = info->sketch->width * info->sketch->depth * sizeof(int) + sizeof(CM_type);
	uint cold_sketch = cold_sketch_->width * cold_sketch_->depth * sizeof(int) + sizeof(CM_type);
	uint cold_heap = cold_heap_->GetSize();
	uint block_index = eps_1_ * eps_1_ * eps_1_ * sizeof(block);
	uint node_count = 0;
	for (int i = 0; i < (int) info->tree_node_count.size(); i++) {
		node_count += info->tree_node_count[i];
	}
	uint node_size = 2 * sizeof(void*) + sizeof(int) + sizeof(RangeTree::Point<int,int>);  // 2 pointers, node_key and sizeof a point
	uint rtree = rtrees.size() * (sizeof(int) + sizeof(void*)) + node_count * node_size;  // trees are extra pointer and tree id + node size
	uint reservoir = reservoir_->getSize();
	if (cdarq) reservoir = 0;

	return tree_sketch + cold_sketch + cold_heap + block_index + rtree + reservoir;
}


int MARQ::GetNrHotCells() {
	int count = 0;
	for (int i = 0; i < eps_1_; i++) {
		for (int j = 0; j < eps_1_; j++) {
			for (int k = 0; k < eps_1_; k++) {
				if (block_index_[i][j][k].nearest_core != -1) {
					count++;
					//std::cout << i << ", " << j << ", " << k << std::endl;
				}
			}
		}
	}
	return count;
}


// Compute block number <x, y, z, ..> from coordinates <x, y, z, ..>
std::vector<int> MARQ::fromCorToBlock(std::vector<int> cor) {
	std::vector<int> result;

	for (int i = 0; i < d_; i++) {
		int step = floor(float(upper_limits_[i])/float(eps_1_));
		int p = ceil(float(cor[i])/float(step));
		result.push_back(p);
	}

	return result;
}


int MARQ::FindHotNeighbour(int x, int y, int z) {
	if (x + 1 > eps_1_ && block_index_[x-1][y][z].nearest_core != -1) {
		return block_index_[x-1][y][z].nearest_core;
	}
	if (x - 1 > 0 && block_index_[x-1][y][z].nearest_core != -1) {
		return block_index_[x-1][y][z].nearest_core;
	}
	if (y + 1 < eps_1_ && block_index_[x][y+1][z].nearest_core != -1) {
		return block_index_[x][y+1][z].nearest_core;
	}
	if (y - 1 > 0 && block_index_[x][y-1][z].nearest_core != -1) {
		return block_index_[x][y-1][z].nearest_core;
	}
	if (z + 1 < eps_1_ && block_index_[x][y][z+1].nearest_core != -1) {
		return block_index_[x][y][z+1].nearest_core;
	}
	if (z - 1 > 0 && block_index_[x][y][z-1].nearest_core != -1) {
		return block_index_[x][y][z-1].nearest_core;
	}
	return -1;	
}


void MARQ::update(std::vector<int> x) {
	// COmpute block
	int firstCor = floor (double(x[0]) / double(ceil(double(upper_limits_[0])/double(eps_1_))));
	int secondCor = floor (double(x[1]) / double (ceil(double(upper_limits_[1])/double(eps_1_))));
	int third = floor (double((uint) SignedIPToUnsigned(x[2]))) / double (ceil(double(UINT_MAX/*eps_1_*/)/double(eps_1_)));  // hard code fix for uint ip
	
	float exactx = (float(x[0]) / float(ceil(float(upper_limits_[0])/float(eps_1_))));
	float exacty = (float(x[1]) / float (ceil(float(upper_limits_[1])/float(eps_1_))));
	float exactz = (float((uint) SignedIPToUnsigned(x[2]))) / float (ceil(float(UINT_MAX/*eps_1_*/)/float(eps_1_)));
	
	// Sample point in blocks
	if (!cdarq) {
		reservoir_->insert(firstCor, secondCor, third, marq_item(exactx, exacty, exactz));
	}

	block &block = block_index_[firstCor][secondCor][third];
	//std::cout << x[0] << ", " << x[1] << ", (" << x[2] << ") " << (uint) SignedIPToUnsigned(x[2]) << " to " << firstCor << ", " << secondCor << ", " << third << std::endl;

	// If block is part of core group
	N_++;  // Increase total number of data points
 	if (block_index_[firstCor][secondCor][third].nearest_core != -1) {
		// Update current block in its range tree, updating all nodes along the way
		rtrees[block.nearest_core].updateSketch(new RangeTree::Point<int,int> (f3(firstCor,secondCor,third), 0)); 
	} else {
		block_index_[firstCor][secondCor][third].count++;  // increase counter of block index
		ccold_++;  // increase cold count

		// Update or insert block into coldheap
		UpdateCold(firstCor, secondCor, third);

		// Cold threshold check
		if ((float) ccold_ > (1 - theta_) * N_) {
			// Get max block from maxheap
			heap_entry b = cold_heap_->extractMax();
			//block_index_[b.x][b.y][b.z].core = true; // Set block to core

			int core_group = FindHotNeighbour(b.x, b.y, b.z);
			if (core_group != -1) { // If there is a core group near this group
				ColdToExistingTree(core_group, b);
			} else {
				ColdToNewTree(b);
			}
			ccold_ -= b.count;

			CM_Update(cold_sketch_, BlockToKey(b.x, b.y, b.z), -1*b.count);  // subtract val from countmin
		}
	}
}


// BUlk insert, assume all entries are cold sketch, have a finalize function which afterwards builds the trees and transfers this to hot
void MARQ::updateBulk(std::vector<int> x) {
	// COmpute block
	int firstCor = floor (double(x[0]) / double(ceil(double(upper_limits_[0])/double(eps_1_))));
	int secondCor = floor (double(x[1]) / double (ceil(double(upper_limits_[1])/double(eps_1_))));
	int third = floor (double((uint) SignedIPToUnsigned(x[2]))) / double (ceil(double(UINT_MAX)/double(eps_1_)));  // hard code fix for uint ip
	
	block &block = block_index_[firstCor][secondCor][third];

	// If block is part of core group
	N_++;  // Increase total number of data points
	block_index_[firstCor][secondCor][third].count++;  // increase counter of block index
	ccold_++;  // increase cold count

	// Update or insert block into coldheap
	UpdateCold(firstCor, secondCor, third);
}


void MARQ::finalizeBulk() {
	// Loop until condition is satisfied
	while ((float) ccold_ >= (1 - theta_) * N_) {
		// Keep getting the new hot item
		heap_entry b = cold_heap_->extractMax();
		// Check if neighbour is already core
		int core_group = FindHotNeighbour(b.x, b.y, b.z);
		// If so, add point ot tree but do not update tree sketch
		if (core_group != -1) { 
			auto tree_points = rtrees[core_group].originalPoints; // get points of the tree
			tree_points.push_back(RangeTree::Point<int,int> (f3(b.x, b.y, b.z), 0));
			rtrees[core_group] = RangeTree::RangeTree<int,int>(tree_points, info, core_group); // rebuild tree
			block_index_[b.x][b.y][b.z].nearest_core = core_group; // set core group to nearest core group
		} else {
			// Otherwise create new tree, but do not update tree sketch
			block_index_[b.x][b.y][b.z].nearest_core = rtrees.size(); // set core group to itself
			info->tree_node_count.push_back(0); // add new tree node count
			// Create new range tree with current block
			RangeTree::Point<int, int> p(f3(b.x, b.y, b.z), 0);
			std::vector<RangeTree::Point<int,int>> points = {p};
			rtrees.push_back(RangeTree::RangeTree<int,int>(points, info, info->tree_count));
			info->tree_count++;
		}
		// Remove counts from cold part
		ccold_ -= b.count;
		CM_Update(cold_sketch_, BlockToKey(b.x, b.y, b.z), -1*b.count);  // subtract val from countmin
	}

	// For every block, if it is hot, update the tree sketch.
	for (int x = 0; x < eps_1_; x++) {
		for (int y = 0; y < eps_1_; y++) {
			for (int z = 0; z < eps_1_; z++) {
				if (block_index_[x][y][z].nearest_core != -1) {
					rtrees[block_index_[x][y][z].nearest_core].updateSketch(new RangeTree::Point<int,int> (f3(x, y, z), 0), block_index_[x][y][z].count); // update tree with current count
				}
			}
		}
	}
}


void MARQ::ColdToExistingTree(int core_group, heap_entry b) {
	//std::cout << "insert tree " << core_group << ", " << b.x << ", " << b.y << ", " << b.z << std::endl;
	auto tree_points = rtrees[core_group].originalPoints; // get points of the tree
	tree_points.push_back(RangeTree::Point<int,int> (f3(b.x, b.y, b.z), 0));
	//delete rtrees[core_group];  // delete existing tree
	rtrees[core_group] = RangeTree::RangeTree<int,int>(tree_points, info, core_group); // rebuild tree
	rtrees[core_group].updateSketch(new RangeTree::Point<int,int> (f3(b.x, b.y, b.z), 0), b.count); // update tree with current count
	block_index_[b.x][b.y][b.z].nearest_core = core_group; // set core group to nearest core group
}

void MARQ::ColdToNewTree(heap_entry b) {
	//std::cout << "new tree " << info->tree_count << ", " << b.x << ", " << b.y << ", " << b.z << std::endl;
	block_index_[b.x][b.y][b.z].nearest_core = rtrees.size(); // set core group to itself
	info->tree_node_count.push_back(0); // add new tree node count
	// Create new range tree with current block
	RangeTree::Point<int, int> p(f3(b.x, b.y, b.z), 0);
	std::vector<RangeTree::Point<int,int>> points = {p};
	rtrees.push_back(RangeTree::RangeTree<int,int>(points, info, info->tree_count));
	rtrees.back().updateSketch(new RangeTree::Point<int,int> (f3(b.x, b.y, b.z), 0), b.count);  // bulk update node
	info->tree_count++;
}

void MARQ::UpdateCold(int firstCor, int secondCor, int third) {
	if (block_index_[firstCor][secondCor][third].count == 1) {
		cold_heap_->insertKey(heap_entry(1, firstCor, secondCor, third));  // insert block into heap
	} else {
		int key = cold_heap_->findIndex(heap_entry(0, firstCor, secondCor, third));
		cold_heap_->increaseKey(key);  // find index of block in heap
	}
	// Update cold sketch
	CM_Update(cold_sketch_, BlockToKey(firstCor, secondCor, third), 1);
}

int MARQ::QuerySurplus(int x1, int y1, int z1, int x2, int y2, int z2) {
	int surplus = 0;
	int old_surplus = 0;
	if (surplus_query_method_ == 0) {
		surplus = CM_PointEst(cold_sketch_, BlockToKey(x1, y1, z1));
	} else if (surplus_query_method_ == 1) {
		for (int x = x1; x <= x2; x++) {
			for (int y = y1; y <= y2; y++) {
				for (int z = z1; z <= z2; z++) {
					surplus += CM_PointEst(cold_sketch_, BlockToKey(x, y, z));
					if (surplus < old_surplus) {
						std::cout << "";
					}
					old_surplus = surplus;
				}
			}
		}
	} else if (surplus_query_method_ == 2) {
		for (int x = x1; x <= x2; x++) {
			for (int y = y1; y <= y2; y++) {
				for (int z = z1; z <= z2; z++) {
					surplus += block_index_[x][y][z].count;
				}
			}
		}
	}
	return surplus;
}


// Exact count 3d
long MARQ::countQuery(std::vector<int>& lower, std::vector<int>& upper) {
	int firstCorLower = ceil (double(lower[0]) / double (ceil(double(upper_limits_[0])/double(eps_1_))));
	int secondCorLower = ceil (double(lower[1]) / double (ceil(double(upper_limits_[1])/double(eps_1_))));
	int thirdCorLower = ceil (double((uint) SignedIPToUnsigned(lower[2]))) / double (ceil(double(UINT_MAX/*eps_1_*/)/double(eps_1_)));

	int firstCorUpper = floor (double(upper[0]) / double (ceil(double(upper_limits_[0])/double(eps_1_))));
	int secondCorUpper = floor (double(upper[1]) / double (ceil(double(upper_limits_[1])/double(eps_1_))));
	int thirdCorUpper = floor (double((uint) SignedIPToUnsigned(upper[2]))) / double (ceil(double(UINT_MAX/*eps_1_*/)/double(eps_1_)));

	float floatLowerX = (float(lower[0]) / float(ceil(float(upper_limits_[0])/float(eps_1_))));
	float floatLowerY = (float(lower[1]) / float (ceil(float(upper_limits_[1])/float(eps_1_))));
	float floatLowerZ = (float((uint) SignedIPToUnsigned(lower[2]))) / float (ceil(float(UINT_MAX/*eps_1_*/)/float(eps_1_)));

	float floatUpperX = (float(upper[0]) / float(ceil(float(upper_limits_[0])/float(eps_1_))));
	float floatUpperY = (float(upper[1]) / float (ceil(float(upper_limits_[1])/float(eps_1_))));
	float floatUpperZ = (float((uint) SignedIPToUnsigned(upper[2]))) / float (ceil(float(UINT_MAX/*eps_1_*/)/float(eps_1_)));

	long result = 0;
	long cold_surplus = 0;
	long sample_surplus = 0;
	// Only query blocks if query ranges cover full blocks5
	if ((floatUpperX - floatLowerX >= 1.0f && floatUpperY - floatLowerY >= 1.0f && floatUpperZ - floatLowerZ >= 1.0f) || cdarq) {
		for (auto tree : rtrees) {
			result += tree.countInRange(f3(firstCorLower,secondCorLower, thirdCorLower), f3(firstCorUpper,secondCorUpper, thirdCorUpper));
		}

		cold_surplus = QuerySurplus(firstCorLower, secondCorLower, thirdCorLower, firstCorUpper, secondCorUpper, thirdCorUpper);
	}
	if (!cdarq) {
		sample_surplus = SampleQuery(lower, upper);
	}

	//std::cout << "rtree result: " << result << ", cold: " << cold_surplus << ", sample: " << sample_surplus << "\n";

	return result + cold_surplus + sample_surplus;
}


int MARQ::SampleQuery(std::vector<int>& lower, std::vector<int>& upper) {
	int exactLowerX = floor (double(lower[0]) / double (ceil(double(upper_limits_[0])/double(eps_1_))));
	int exactLowerY = floor (double(lower[1]) / double (ceil(double(upper_limits_[1])/double(eps_1_))));
	int exactLowerZ = floor (double((uint) SignedIPToUnsigned(lower[2]))) / double (ceil(double(UINT_MAX/*eps_1_*/)/double(eps_1_)));

	int exactUpperX = floor (double(upper[0]) / double (ceil(double(upper_limits_[0])/double(eps_1_))));
	int exactUpperY = floor (double(upper[1]) / double (ceil(double(upper_limits_[1])/double(eps_1_))));
	int exactUpperZ = floor (double((uint) SignedIPToUnsigned(upper[2]))) / double (ceil(double(UINT_MAX/*eps_1_*/)/double(eps_1_)));

	float floatLowerX = (float(lower[0]) / float(ceil(float(upper_limits_[0])/float(eps_1_))));
	float floatLowerY = (float(lower[1]) / float (ceil(float(upper_limits_[1])/float(eps_1_))));
	float floatLowerZ = (float((uint) SignedIPToUnsigned(lower[2]))) / float (ceil(float(UINT_MAX/*eps_1_*/)/float(eps_1_)));

	float floatUpperX = (float(upper[0]) / float(ceil(float(upper_limits_[0])/float(eps_1_))));
	float floatUpperY = (float(upper[1]) / float (ceil(float(upper_limits_[1])/float(eps_1_))));
	float floatUpperZ = (float((uint) SignedIPToUnsigned(upper[2]))) / float (ceil(float(UINT_MAX/*eps_1_*/)/float(eps_1_)));


	std::vector<std::vector<std::vector<bool>>> visited(eps_1_, std::vector<std::vector<bool>>(eps_1_, std::vector<bool>(eps_1_, false)));

	int estimate = 0;

	// unrounded query range
	std::vector<float> lower_range = { floatLowerX, floatLowerY, floatLowerZ };
	std::vector<float> upper_range = { floatUpperX, floatUpperY, floatUpperZ };

	// x margin
	if (floatLowerX < (float) exactLowerX || floatUpperX > (float) exactUpperX) {
		for (int y = floor(floatLowerY); y <= floor(floatUpperY); y++) {
			for (int z = floor(floatLowerZ); z <= floor(floatUpperZ); z++) {
				// if block not visited
				if (floatLowerX < (float) exactLowerX) {
					if (!visited[exactLowerX][y][z]) {
						// mark block as visited
						visited[exactLowerX][y][z] = true;
						// Get intersection of samples and exact count
						std::vector<marq_item> unioned_samples = reservoir_->getCombinedSamples(exactLowerX, y, z);
						if (unioned_samples.size() == 0) {
							continue;
						}
						int item_in_sample_frequency = reservoir_->getFilteredUnionFrequency(unioned_samples, lower_range, upper_range);
						int exact_frequency = block_index_[exactLowerX][y][z].count;
						// scale up estimate
						estimate += (int) (item_in_sample_frequency / (float) unioned_samples.size()) * exact_frequency;
					}
				} else if (floatUpperX > (float) exactUpperX) {
					if (!visited[exactUpperX][y][z]) {
						// mark block as visited
						visited[exactUpperX][y][z] = true;
						// Get intersection of samples and exact count
						std::vector<marq_item> unioned_samples = reservoir_->getCombinedSamples(exactUpperX, y, z);
						if (unioned_samples.size() == 0) {
							continue;
						}
						int item_in_sample_frequency = reservoir_->getFilteredUnionFrequency(unioned_samples, lower_range, upper_range);
						int exact_frequency = block_index_[exactUpperX][y][z].count;
						// scale up estimate
						estimate += (int) (item_in_sample_frequency / (float) unioned_samples.size()) * exact_frequency;
					}
				}
			}
		}
	}

	// y margin
	if (floatLowerY < (float) exactLowerY || floatUpperY > (float) exactUpperY) {
		for (int x = floor(floatLowerX); x <= floor(floatUpperX); x++) {
			for (int z = floor(floatLowerZ); z <= floor(floatUpperZ); z++) {
				// if block not visited
				if (floatLowerY < (float) exactLowerY) {
					if (!visited[x][exactLowerY][z]) {
						// mark block as visited
						visited[x][exactLowerY][z] = true;
						// Get intersection of samples and exact count
						std::vector<marq_item> unioned_samples = reservoir_->getCombinedSamples(x, exactLowerY, z);
						if (unioned_samples.size() == 0) {
							continue;
						}
						int item_in_sample_frequency = reservoir_->getFilteredUnionFrequency(unioned_samples, lower_range, upper_range);
						int exact_frequency = block_index_[x][exactLowerY][z].count;
						// scale up estimate
						estimate += (int) (item_in_sample_frequency / (float) unioned_samples.size()) * exact_frequency;
					}
				} else if (floatUpperY > (float) exactUpperY) {
					if (!visited[x][exactUpperY][z]) {
						// mark block as visited
						visited[x][exactUpperY][z] = true;
						// Get intersection of samples and exact count
						std::vector<marq_item> unioned_samples = reservoir_->getCombinedSamples(x, exactUpperY, z);
						if (unioned_samples.size() == 0) {
							continue;
						}
						int item_in_sample_frequency = reservoir_->getFilteredUnionFrequency(unioned_samples, lower_range, upper_range);
						int exact_frequency = block_index_[x][exactUpperY][z].count;
						// scale up estimate
						estimate += (int) (item_in_sample_frequency / (float) unioned_samples.size()) * exact_frequency;
					}
				}
			}
		}
	}

	// z margin
	if (floatLowerZ < (float) exactLowerZ || floatUpperZ > (float) exactUpperZ) {
		for (int x = floor(floatLowerX); x <= floor(floatUpperX); x++) {
			for (int y = floor(floatLowerY); y <= floor(floatUpperY); y++) {
				// if block not visited
				if (floatLowerZ < (float) exactLowerZ) {
					if (!visited[x][y][exactLowerZ]) {
						// mark block as visited
						visited[x][y][exactLowerZ] = true;
						// Get intersection of samples and exact count
						std::vector<marq_item> unioned_samples = reservoir_->getCombinedSamples(x, y, exactLowerZ);
						if (unioned_samples.size() == 0) {
							continue;
						}
						int item_in_sample_frequency = reservoir_->getFilteredUnionFrequency(unioned_samples, lower_range, upper_range);
						int exact_frequency = block_index_[x][y][exactLowerZ].count;
						// scale up estimate
						estimate += (int) (item_in_sample_frequency / (float) unioned_samples.size()) * exact_frequency;
					}
				} else if (floatUpperZ > (float) exactUpperZ) {
					if (!visited[x][y][exactUpperZ]) {
						// mark block as visited
						visited[x][y][exactUpperZ] = true;
						// Get intersection of samples and exact count
						std::vector<marq_item> unioned_samples = reservoir_->getCombinedSamples(x, y, exactUpperZ);
						if (unioned_samples.size() == 0) {
							continue;
						}
						int item_in_sample_frequency = reservoir_->getFilteredUnionFrequency(unioned_samples, lower_range, upper_range);
						int exact_frequency = block_index_[x][y][exactUpperZ].count;
						// scale up estimate
						estimate += (int) (item_in_sample_frequency / (float) unioned_samples.size()) * exact_frequency;
					}
				}
			}
		}
	}


	return estimate;
}