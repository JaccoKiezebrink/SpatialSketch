#ifndef _ELASTIC_SKETCH_H_
#define _ELASTIC_SKETCH_H_

#include "HeavyPart.h"
#include "LightPart.h"
#include "param.h"



template<int bucket_num, int tot_memory_in_bytes>
class ElasticSketch
{
    static constexpr int heavy_mem = bucket_num * COUNTER_PER_BUCKET * 8;
    static constexpr int light_mem = tot_memory_in_bytes - heavy_mem;
    
    BOBHash32 *bobhash = NULL;
    HeavyPart<bucket_num> heavy_part;
    LightPart<light_mem> light_part;

public:
    ElasticSketch(){
        std::random_device rd;
        bobhash = new BOBHash32(rd() % MAX_PRIME32);

        // hashab_ = new uint32_t[3];
        // std::srand(7);
        // //std:: cout << std::rand() << ", " << std::rand() << ", " << std::rand() << std::endl;
        // //hashab_[i][0] = int(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
        // //hashab_[i][1] = int(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
        // //hashab_[i][0] = int(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
        // hashab_[0] = std::rand();//long(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
        // hashab_[1] = std::rand();//long(float(std::rand())*float(92821)/float(RAND_MAX) + 1);
        // // random uint32_t:

        // // generate random large prime number
        // long prime = 10000 + ( std::rand() % ( RAND_MAX - 10000 + 1 ) );
        // int count = 0;
        // while (!ElasticSketch<bucket_num, tot_memory_in_bytes>::checkIfPrime(prime) && count < 1000) {
        //     prime = 10000 + ( std::rand() % ( RAND_MAX - 10000 + 1 ) );
        // // }
        
        // hashab_[2] = prime;
}
    ~ElasticSketch(){
        delete bobhash;
    }
    void clear();

    void insert(uint8_t *key, int f = 1);
    void insert(uint8_t *key, int &f, uint* hashes);
    void quick_insert(uint8_t *key, int f = 1);

    int query(uint8_t *key);
    int query_compressed_part(uint8_t *key, uint8_t *compress_part, int compress_counter_num);

    int get_compress_width(int ratio) { return light_part.get_compress_width(ratio);}
    void compress(int ratio, uint8_t *dst) {    light_part.compress(ratio, dst); }

    int get_bucket_num() { return heavy_part.get_bucket_num(); }
    double get_bandwidth(int compress_ratio) ;

    void get_heavy_hitters(int threshold, vector<pair<string, int>> & results);
    int get_cardinality();
    double get_entropy();
    void get_distribution(vector<double> &dist);

    void *operator new(size_t sz);
    void operator delete(void *p);

    int GetItemHashes(uint8_t *key, uint32_t* hashes); // Compute hashes
    bool checkIfPrime(long prime);
private:
        uint32_t *hashab_;
};

template<int bucket_num, int tot_memory_in_bytes>
bool ElasticSketch<bucket_num,tot_memory_in_bytes>::checkIfPrime(long prime) {
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

template<int bucket_num, int tot_memory_in_bytes>
void ElasticSketch<bucket_num, tot_memory_in_bytes>::clear()
{
    heavy_part.clear();
    light_part.clear();
   
}

// Compute hashes
template<int bucket_num, int tot_memory_in_bytes>
int ElasticSketch<bucket_num, tot_memory_in_bytes>::GetItemHashes(uint8_t *key, uint32_t* hashes) {
    //memset(hashes, 0, sizeof(uint) * depth_);
    uint32_t fp = *((uint32_t*)key);
    hashes[0] = CalculateBucketPos(fp) % bucket_num; // Heavy bucket position.
    
    //hashes[1] = (uint32_t) ((hashab_[0] * *key + hashab_[1]) % hashab_[2] % light_mem + light_mem) % light_mem;
    uint32_t hash_val = (uint32_t)bobhash->run((const char*)key, KEY_LENGTH_4) % light_mem;
    hashes[1] = hash_val;
    return 2;
}

// Pass hashes to avoid recalculation
template<int bucket_num, int tot_memory_in_bytes>
inline void ElasticSketch<bucket_num, tot_memory_in_bytes>::insert(uint8_t *key, int &f, uint32_t *hashes)
{
    uint8_t swap_key[KEY_LENGTH_4];
    uint32_t swap_val = 0;
    int result = heavy_part.insert_precomp_hashes(key, swap_key, swap_val, hashes[0], f);
    switch(result)
    {
        case 0: return;
        case 1:{
            if(HIGHEST_BIT_IS_1(swap_val))
                light_part.insert(swap_key, GetCounterVal(swap_val), hashes[1]);
            else
                light_part.swap_insert(swap_key, swap_val, hashes[1]);
            return;
        }
        case 2: {
            light_part.insert(key, f, hashes[1]);
            return;
            }
        default:
            printf("error return value !\n");
            exit(1);
    }
}

template<int bucket_num, int tot_memory_in_bytes>
inline void ElasticSketch<bucket_num, tot_memory_in_bytes>::insert(uint8_t *key, int f)
{
    uint8_t swap_key[KEY_LENGTH_4];
    uint32_t swap_val = 0;
    int result = heavy_part.insert(key, swap_key, swap_val, f);
    switch(result)
    {
        case 0: return;
        case 1:{
            if(HIGHEST_BIT_IS_1(swap_val))
                light_part.insert(swap_key, GetCounterVal(swap_val));
            else
                light_part.swap_insert(swap_key, swap_val);
            return;
        }
        case 2: light_part.insert(key, 1);  return;
        default:
            printf("error return value !\n");
            exit(1);
    }
}

template<int bucket_num, int tot_memory_in_bytes>
void ElasticSketch<bucket_num, tot_memory_in_bytes>::quick_insert(uint8_t *key, int f)
{
    heavy_part.quick_insert(key, f);
}

template<int bucket_num, int tot_memory_in_bytes>
inline int ElasticSketch<bucket_num, tot_memory_in_bytes>::query(uint8_t *key)
{
    uint32_t heavy_result = heavy_part.query(key);
    if(heavy_result == 0 || HIGHEST_BIT_IS_1(heavy_result))
    {
        int light_result = light_part.query(key);
        return (int)GetCounterVal(heavy_result) + light_result;
    }
    return heavy_result;
}

template<int bucket_num, int tot_memory_in_bytes>
int ElasticSketch<bucket_num, tot_memory_in_bytes>::query_compressed_part(uint8_t *key, uint8_t *compress_part, int compress_counter_num)
{
    uint32_t heavy_result = heavy_part.query(key);
    if(heavy_result == 0 || HIGHEST_BIT_IS_1(heavy_result))
    {
        int light_result = light_part.query_compressed_part(key, compress_part, compress_counter_num);
        return (int)GetCounterVal(heavy_result) + light_result;
    }
    return heavy_result;
}

template<int bucket_num, int tot_memory_in_bytes>
double ElasticSketch<bucket_num, tot_memory_in_bytes>::get_bandwidth(int compress_ratio) 
{
    int result = heavy_part.get_memory_usage();
    result += get_compress_width(compress_ratio) * sizeof(uint8_t);
    return result * 1.0 / 1024 / 1024;
}

template<int bucket_num, int tot_memory_in_bytes>
void ElasticSketch<bucket_num, tot_memory_in_bytes>::get_heavy_hitters(int threshold, vector<pair<string, int>> & results)
{
    for (int i = 0; i < bucket_num; ++i) 
        for (int j = 0; j < MAX_VALID_COUNTER; ++j) 
        {
            uint32_t key = heavy_part.buckets[i].key[j];
            int val = query((uint8_t *)&key);
            if (val >= threshold) {
                results.push_back(make_pair(string((const char*)&key, 4), val));
            }
        }
}

template<int bucket_num, int tot_memory_in_bytes>
int ElasticSketch<bucket_num, tot_memory_in_bytes>::get_cardinality()
{
    int card = light_part.get_cardinality();
    for(int i = 0; i < bucket_num; ++i)
        for(int j = 0; j < MAX_VALID_COUNTER; ++j)
        {
            uint8_t key[KEY_LENGTH_4];
            *(uint32_t*)key = heavy_part.buckets[i].key[j];
            int val = heavy_part.buckets[i].val[j];
            int ex_val = light_part.query(key);

            if(HIGHEST_BIT_IS_1(val) && ex_val)
            {
                val += ex_val;
                card--;
            }
            if(GetCounterVal(val))
                card++;
        }
    return card;
}

template<int bucket_num, int tot_memory_in_bytes>
double ElasticSketch<bucket_num, tot_memory_in_bytes>::get_entropy()
{
    int tot = 0;
    double entr = 0;

    light_part.get_entropy(tot, entr);

    for(int i = 0; i < bucket_num; ++i)
        for(int j = 0; j < MAX_VALID_COUNTER; ++j)
        {
            uint8_t key[KEY_LENGTH_4];
            *(uint32_t*)key = heavy_part.buckets[i].key[j];
            int val = heavy_part.buckets[i].val[j];

            int ex_val = light_part.query(key);

            if(HIGHEST_BIT_IS_1(val) && ex_val)
            {
                val += ex_val;

                tot -= ex_val;

                entr -= ex_val * log2(ex_val);
            }
            val = GetCounterVal(val);
            if(val)
            {
                tot += val;
                entr += val * log2(val);
            }
        }
    return -entr / tot + log2(tot);
}

template<int bucket_num, int tot_memory_in_bytes>
void ElasticSketch<bucket_num, tot_memory_in_bytes>::get_distribution(vector<double> &dist)
{
    light_part.get_distribution(dist);

    for(int i = 0; i < bucket_num; ++i)
        for(int j = 0; j < MAX_VALID_COUNTER; ++j)
        {
            uint8_t key[KEY_LENGTH_4];
            *(uint32_t*)key = heavy_part.buckets[i].key[j];
            int val = heavy_part.buckets[i].val[j];

            int ex_val = light_part.query(key);

            if(HIGHEST_BIT_IS_1(val) && ex_val != 0)
            {
                val += ex_val;
                dist[ex_val]--;
            }
            val = GetCounterVal(val);
            if(val)
            {
                if(val + 1 > dist.size())
                    dist.resize(val + 1);
                dist[val]++;
            }
        }
}

template<int bucket_num, int tot_memory_in_bytes>
void* ElasticSketch<bucket_num, tot_memory_in_bytes>::operator new(size_t sz)
{
    constexpr uint32_t alignment = 64;
    size_t alloc_size = (2 * alignment + sz) / alignment * alignment;
    void *ptr = ::operator new(alloc_size);
    void *old_ptr = ptr;
    void *new_ptr = ((char*)std::align(alignment, sz, ptr, alloc_size) + alignment);
    ((void **)new_ptr)[-1] = old_ptr;

    return new_ptr;
}

template<int bucket_num, int tot_memory_in_bytes>
void ElasticSketch<bucket_num, tot_memory_in_bytes>::operator delete(void *p)
{
    ::operator delete(((void**)p)[-1]);
}


#endif // _ELASTIC_SKETCH_H_
