#ifndef UTILS_H_
#define UTILS_H_

#include <vector>
#include <string>
#include <memory>
#include <limits.h>

// Utility structs, functions, etc. to be used in various places

typedef struct rect {
    float x1, y1, x2, y2;
} rect;

// Query have to be integers
typedef struct range {
    int x1, y1, x2, y2;

    range() {}
    range(int x1, int y1, int x2, int y2) : x1(x1), y1(y1), x2(x2), y2(y2) {}
} range;


typedef struct shape_info {
    int grid_size = -1;
    int selection_size = -1;
    std::string shape = "";
    int max_x_offset = 0; 
    int max_y_offset = 0;
    int hole_vertex_state = 0; // 0: normal vertex, 1: hole vertex, 2: line vertex connecting main polygon to hole
    std::vector<std::pair<float, float>> vertices, hole_vertices;
    std::vector<std::pair<int, int>> coordinates, offset_coordinates;
} shape_info;

// Convert rectangle to range
range RectToRange(rect r);

// Parsing of vertices or coordinates
std::vector<std::pair<int, int>> ParseIntPairs(std::ifstream &file);
std::vector<std::pair<float, float>> parse_float_pairs(std::ifstream &file);

bool ParseShapeFile(std::string file_name, shape_info &shape, bool subtractive_querying=false);

// Check range is within NxN grid
bool RangeBoundsCheck(range &rang, int N);

// Given nr of bytes, return scaled value to closest order unit (KB, MB, GB, etc.)
std::pair<float, std::string> ScaleMemory(uint bytes);
std::pair<float, std::string> ScaleMemory(float bytes);

// Compare functions for sorting
struct CompFirstFloatPair {
    bool operator()(const float &a, const std::pair<float,float> &b) const { return a < b.first; }
    bool operator()(const std::pair<float,float> &a, const float &b) const { return a.first < b; }
};

struct CompSecondFloatPair {
    bool operator()(const float &a, const std::pair<float,float> &b) const { return a < b.second; }
    bool operator()(const std::pair<float,float> &a, const float &b) const { return a.second < b; }
};

struct CompFirstFloatTuple {
    bool operator()(const float &a, const std::tuple<float, float, int> &b) const { return a < std::get<0>(b); }
    bool operator()(const std::tuple<float, float, int>&a, const float &b) const { return std::get<0>(a) < b; }
};

struct CompSecondFloatTuple {
    bool operator()(const float &a, const std::tuple<float, float, int> &b) const { return a < std::get<1>(b); }
    bool operator()(const std::tuple<float, float, int> &a, const float &b) const { return std::get<1>(a) < b; }
};

bool CompBySecondInt(const std::pair<int,int> &a, const std::pair<int,int> &b);
bool CompByFirstFloat(const std::pair<float,float> &a, const std::pair<float,float> &b);
bool CompByFirstTupleFloat(const std::tuple<float, float, int> &a, const std::tuple<float, float, int> &b);
bool CompBySecondFloat(const std::pair<float,float> &a, const std::pair<float,float> &b);
bool CompBySecondTupleFloat(const std::tuple<float, float, int> &a, const std::tuple<float, float, int> &b);

int SignedIPToUnsigned(uint ip);

uint UnsignedIPToSigned(int ip);

bool isPrime(long prime);

#endif  // UTILS_H_