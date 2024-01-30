#include "Utils.h"

#include <iostream>
#include <fstream>
#include <cmath>

//#define DEBUG  // shape file parsing debug printing

// Rectangles from the partitioning lay on an offset of 0.5, these can be converted to ranges of integers
range RectToRange(rect r) {
    range range;
    range.x1 = (int) std::ceil(r.x1);
    range.y1 = (int) std::ceil(r.y1);
    range.x2 = (int) std::floor(r.x2);
    range.y2 = (int) std::floor(r.y2);
    return range;
}


// Given open file parse pairs of floats delimited by , and end line \n
std::vector<std::pair<int, int>> ParseIntPairs(std::ifstream &file) {
    std::vector<std::pair<int, int>> pairs;

    std::string line;
    std::string delim = ",";
    while (std::getline(file, line)) {
        auto pos = line.find(",");
        if (pos == std::string::npos) {
            break;
        }
        std::string x = line.substr(0, pos);
        line.erase(0, line.find(delim) + delim.length());
        std::string y = line;       
        pairs.push_back(std::make_pair(std::stoi(x), std::stof(y)));
    }

    return pairs;
}

std::vector<std::pair<float, float>> parse_float_pairs(std::ifstream &file) {
    std::vector<std::pair<float, float>> pairs;

    std::string line;
    std::string delim = ",";
    while (std::getline(file, line)) {
        // Skip hole related things
        if (line.find("hole") == std::string::npos) {
            continue;
        }
        auto pos = line.find(",");
        if (pos == std::string::npos) {
            break;
        }
        std::string x = line.substr(0, pos);
        line.erase(0, line.find(delim) + delim.length());
        std::string y = line;       
        pairs.push_back(std::make_pair(std::stof(x), std::stof(y)));
    }

    return pairs;
}


bool ParseShapeFile(std::string file_name, shape_info &shape_info, bool subtractive_querying) {
    std::ifstream file;
    file.open(file_name);
    if (!file.is_open()) {
        std::cout << "Error opening file: " << file_name << std::endl;
        return false;
    }

    std::vector<std::pair<float, float>> vertices, hole_vertices;
    std::vector<std::pair<int, int>> coordinates, offset_coordinates;

    std::string line;
    while (std::getline(file, line)) {
        if (line.find("grid size:") != std::string::npos) {
            shape_info.grid_size = std::stoi(line.substr(line.find(":") + 1));
#ifdef DEBUG
            std::cout << "grid size: " << shape.grid_size << std::endl;
#endif
        } else if (line.find("selection size:") != std::string::npos) {
            shape_info.selection_size = std::stoi(line.substr(line.find(":") + 1));
#ifdef DEBUG
            std::cout << "selection_size: " << shape.selection_size << std::endl;
#endif
        } else  if (line.find("shape:") != std::string::npos) {
            shape_info.shape = line.substr(line.find(":") + 1);
            // remove endline
            shape_info.shape.erase(std::remove(shape_info.shape.begin(), shape_info.shape.end(), '\r'), shape_info.shape.end());

#ifdef DEBUG
            std::cout << "shape:" << shape << std::endl;
#endif
        } else  if (line.find("shape_type:") != std::string::npos) {
            std::string shape_type = line.substr(line.find(":") + 1);
            // remove endline
            shape_type.erase(std::remove(shape_type.begin(), shape_type.end(), '\r'), shape_type.end());

#ifdef DEBUG
            std::cout << "shape_type:" << shape_type << std::endl;
#endif
        } else  if (line.find("max_x_offset:") != std::string::npos) {
            shape_info.max_x_offset = std::max(std::stoi(line.substr(line.find(":") + 1)), shape_info.max_x_offset);
#ifdef DEBUG
            std::cout << "max_x_offset:" << max_x_offset << std::endl;
#endif
        } else  if (line.find("max_y_offset:") != std::string::npos) {
            std::string shape_type = line.substr(line.find(":") + 1);
            shape_info.max_y_offset = std::max(std::stoi(line.substr(line.find(":") + 1)), shape_info.max_y_offset);
#ifdef DEBUG
            std::cout << "max_y_offset:" << max_y_offset << std::endl;
#endif
        } else if (line.find("vertices") != std::string::npos) {                
            std::string line;
            std::string delim = ",";
            while (std::getline(file, line)) {
                // Check for hole/line start/end
                if (line.find("hole start") != std::string::npos) {
                    // Indicate following vertices are part of hole
                    shape_info.hole_vertex_state = 1;
                    continue;
                } else if (line.find("hole end") != std::string::npos) {
                    // Indicate following vertices are not part of hole anymore, therefore will be line vertex
                    shape_info.hole_vertex_state = 2;
                    continue;
                } else if (line.find("line start") != std::string::npos) {
                    // Indicate following vertices are part of line connecting main polygon to hole
                    shape_info.hole_vertex_state = 2;
                    continue;
                } else if (line.find("line end") != std::string::npos) {
                    // End of line connecting main polygon to hole
                    shape_info.hole_vertex_state = 0;
                    continue;
                }

                auto pos = line.find(",");
                if (pos == std::string::npos) {
                    break;
                }
                
                std::string x = line.substr(0, pos);
                line.erase(0, line.find(delim) + delim.length());
                std::string y = line;       

                // If the current vertex is part of a hole AND we want to do something with this (subtractive querying), then push it to seperate vector
                if (shape_info.hole_vertex_state == 1 && subtractive_querying) {
                    hole_vertices.push_back(std::make_pair(std::stof(x), std::stof(y)));
                } else if (shape_info.hole_vertex_state == 2 && subtractive_querying) { // if line vertex, ignore
                    continue;
                } else { // otherwise keep all vertices together in the same vector
                    vertices.push_back(std::make_pair(std::stof(x), std::stof(y)));
                }
            }
            // hole vertices have to be reversed, as their order is clockwise, but now has to be counter-clockwise 
            std::reverse(hole_vertices.begin(), hole_vertices.end());

#ifdef DEBUG
            std::cout << "vertices " << vertices.size() << std::endl;
#endif
        } else if (line.find("coordinates") != std::string::npos) {
            coordinates = ParseIntPairs(file);
#ifdef DEBUG
            std::cout << "coordinates " << coordinates.size() << std::endl;
#endif
        } else {
            std::cout << "Unexpected line in file: " << line << std::endl;
            return false;
        }
    }

    shape_info.vertices = vertices;
    shape_info.hole_vertices = hole_vertices;
    shape_info.coordinates = coordinates;

    file.close();

    return true;
}


bool RangeBoundsCheck(range &rang, int N) {
    bool correct = true;
    if (rang.x1 < 0 || rang.x2 >= N || rang.y1 < 0 || rang.y2 >= N) {
        //std::cout << "Error: query range (" << rang.x1 << ", " << rang.y1 << ", " << rang.x2 << ", " << rang.y2 << ") out of bounds, ";
        if (rang.x1 < 0) {
            rang.x1 = 0;
            correct = false;
        } 
        if (rang.x2 >= N) {
            rang.x2 = N - 1;
            if (rang.x2 < rang.x1) {
                rang.x1 = rang.x2;
            }
            correct = false;
        }
        if (rang.y1 < 0) {
            rang.y1 = 0;
            correct = false;
        }
        if (rang.y2 >= N) {
            rang.y2 = N - 1;
            if (rang.y2 < rang.y1) {
                rang.y1 = rang.y2;
            }
            correct = false;
        }
        //std::cout << "corrected to (" << rang.x1 << ", " << rang.y1 << ", " << rang.x2 << ", " << rang.y2 << ") out of bounds, " << std::endl;
    }
    return correct;
}


std::pair<float, std::string> ScaleMemory(uint bytes) {
    std::string unit = "B";
    double mem = (double) bytes;
    if (mem > 1024) {
        mem /= 1024;
        unit = "KB";
    }
    if (mem > 1024) {
        mem /= 1024;
        unit = "MB";
    }
    if (mem > 1024) {
        mem /= 1024;
        unit = "GB";
    }
    return std::make_pair((float) mem, unit);
}

std::pair<float, std::string> ScaleMemory(float bytes) {
    std::string unit = "B";
    if (bytes > 1024) {
        bytes /= 1024;
        unit = "KB";
    }
    if (bytes > 1024) {
        bytes /= 1024;
        unit = "MB";
    }
    if (bytes > 1024) {
        bytes /= 1024;
        unit = "GB";
    }
    return std::make_pair(bytes, unit);
}


// Function to sort pairs by second attribute
bool CompBySecondInt(const std::pair<int,int> &a, const std::pair<int,int> &b) {
    return (a.second < b.second);
}

bool CompByFirstFloat(const std::pair<float,float> &a, const std::pair<float,float> &b) {
    return (a.first < b.first);
}

bool CompBySecondFloat(const std::pair<float,float> &a, const std::pair<float,float> &b) {
    return (a.second < b.second);
}

bool CompByFirstTupleFloat(const std::tuple<float, float, int> &a, const std::tuple<float, float, int> &b) {
    return (std::get<0>(a) < std::get<0>(b));
}

bool CompBySecondTupleFloat(const std::tuple<float, float, int> &a, const std::tuple<float, float, int> &b) {
    return (std::get<1>(a) < std::get<1>(b));
}

int SignedIPToUnsigned(uint ip) {
    return int(long(ip) - long(INT_MIN));
}

uint UnsignedIPToSigned(int ip) {
    return uint(long(ip) - long(INT_MIN));
}

