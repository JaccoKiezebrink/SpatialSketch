#include "Partitioner.h"

#include <algorithm>
#include <iostream>
#include <cstring>
#include <unordered_map>
#include <cmath>

//#define THROW_EXCEPTIONS  // Define if exceptions should be thrown or simply print statements set
#define PRINT_EXCEPTIONS

// Vertex comparison function for sorting on x and y coordinate respectively
bool CompVertexOnX(const vertex &a, const vertex &b) {
    if (a.x == b.x && a.y == b.y) {
        return a.id < b.id;
    } else if (a.x == b.x) {
        return a.y < b.y;
    }
    return (a.x < b.x);
}

bool CompVertexOnY(const vertex &a, const vertex &b) {
    if (a.y == b.y && a.x == b.x) {
        return a.id < b.id;
    } else if (a.y == b.y) {
        return a.x < b.x;
    }
    return (a.y < b.y);
}


// Constructor takes vertices as input sorted in anti-clockwise order
Partitioner::Partitioner(std::vector<std::pair<float, float>> vertices) {
    orig_vertices_ = vertices;
    N_ = (int) orig_vertices_.size();
    concave_count_ = 0;
    partition_success_ = true; // assume correctness from the start, but set to false if problems occur

    // SKip any computation or memory reservation if there are no vertices
    if (N_ != 0) {
        vertices_.reserve(orig_vertices_.size());
        if (!ComputeVertexType()) {
            partition_success_ = false;
#ifdef PRINT_EXCEPTIONS
            std::cout << "ComputeVertexType failed" << std::endl;
#endif
        }
    }
}

std::vector<rect> Partitioner::ComputePartitioning() {
    std::vector<rect> result = {};

    if (N_ == 0) {
        return result;
    }

    // Compute chords and return smallest of horizontal or vertical set
    if (ComputeMinChords()) {
        //std::cout << "Min chords computed successfully" << std::endl;
    } else {
#ifdef PRINT_EXCEPTIONS
        std::cout << "Min chords computation failed" << std::endl;
#endif
        partition_success_ = false;
        return result;
    }

    std::vector<vertex> reduced_vertices = ReduceVertices(vertices_);
    return CornerIdentification(reduced_vertices);
}

Direction Partitioner::DetermineDirection(std::pair<float, float> v, std::pair<float, float> w) {
    std::pair<float, float> diff;
    diff.first = std::roundf(w.first - v.first);
    diff.second = std::roundf(w.second - v.second);
    
    // Note the 0.5 implies difference between coordinates should be atleast 0.5 manhatten distance
    if (diff.first == 0 && diff.second >= 0.5) {
        return UP;
    } else if (diff.first == 0 && diff.second <= -0.5) {
        return DOWN;
    } else if (diff.first >= 0.5 && diff.second == 0) {
        return RIGHT;
    } else if (diff.first <= -0.5 && diff.second == 0) {
        return LEFT;
    }

    partition_success_ = false;
    std::string except = "!Error: (" + std::to_string(diff.first) + ", " + std::to_string(diff.second) 
                        + ") is not a valid direction\nCoordinates: (" + std::to_string(v.first) + ", " 
                        + std::to_string(v.second) + "), (" + std::to_string(w.first) + ", " + std::to_string(w.second) + ")";
#ifdef THROW_EXCEPTIONS
    throw std::invalid_argument(except);
#endif
#ifdef PRINT_EXCEPTIONS
    std::cout << except << std::endl;
#endif
    return Direction::DIR_UNKNOWN;
}

Direction Partitioner::DetermineDirection(vertex v, vertex w) {
    return DetermineDirection(std::make_pair(v.x, v.y), std::make_pair(w.x, w.y));
}


inline std::string VertexToString(std::pair<float, float> v) {
    return std::to_string(v.first) + "," + std::to_string(v.second);
}


bool Partitioner::ComputeVertexType() {
    int id_count = 0;
    //std::unordered_map<std::string, int> id_map = std::unordered_map<std::string, int>();
    //id_map.reserve(orig_vertices_.size());

    for (size_t i = 0 ; i < orig_vertices_.size(); i++) {
        std::pair<float, float> v1, v2, v3, diff1, diff2;
        Direction dir1, dir2;  // direction of v1-v2 and v2-v3
        int id;

        // obtain three consecutive vertices
        v1 = orig_vertices_[i];
        v2 = orig_vertices_[(i+1)%orig_vertices_.size()];
        v3 = orig_vertices_[(i+2)%orig_vertices_.size()];  

        vertex_map_.insert(std::make_pair(v2, v3));

        /*auto search = id_map.find(VertexToString(v2));
        if (search == id_map.end()) {
            id = id_count;
            id_count++;
            id_map.insert({VertexToString(v2), id});
        } else {
            id = search->second;
        }*/
        id = id_count;
        id_count++;

        dir1 = DetermineDirection(v1, v2);
        dir2 = DetermineDirection(v2, v3);     

        if (dir1 == dir2) {
            vertices_.push_back(vertex(id, v2.first, v2.second, STRAIGHT, dir1));
        } else if ((dir1 == RIGHT && dir2 == DOWN) || (dir1 == LEFT && dir2 == UP) || (dir1 == UP && dir2 == RIGHT) || (dir1 == DOWN && dir2 == LEFT)) {
            vertices_.push_back(vertex(id, v2.first, v2.second, CONCAVE, dir1));
            concave_count_++;
        } else if ((dir1 == RIGHT && dir2 == UP) || (dir1 == LEFT && dir2 == DOWN) || (dir1 == UP && dir2 == LEFT) || (dir1 == DOWN && dir2 == RIGHT)) {
            vertices_.push_back(vertex(id, v2.first, v2.second, CONVEX, dir1));
        } else {
            std::cout << "\n!Error: (" << dir1 << ", " << dir2 << ") by vertices " << (i)%orig_vertices_.size() << ", " << (i+1)%orig_vertices_.size() << ", " << (i+2)%orig_vertices_.size() << " does not create valid angle" << std::endl;
            std::cout << "Coordinates: (" << v1.first << ", " << v1.second << "), (" << v2.first << ", " << v2.second << "), (" << v3.first << ", " << v3.second << ")" << std::endl;
/*#ifdef THROW_EXCEPTIONS
        throw std::runtime_error("ComputeVertexType failes");
#endif*/
            return false;
        }
    }

    return true;
}


bool Partitioner::AreNeighbors(vertex v, vertex w) {
    return (((v.id + 1) % N_ == w.id || (w.id + 1) % N_ == v.id));
    //return std::abs(v.id - w.id) % N_ == 1;
}


bool Partitioner::ComputeMinChords() {
    std::vector<std::pair<vertex, vertex>> hor_chords, ver_chords;
    std::vector<std::pair<int, int>> hor_chords_id, ver_chords_id;
    hor_chords.reserve(concave_count_);
    ver_chords.reserve(concave_count_);

    Direction dir_prev, dir_cur, dir_next;

    // Obtain vertical chords
    std::vector<vertex> x_sorted = vertices_;
    std::sort(x_sorted.begin(), x_sorted.end(), CompVertexOnX);

    for (size_t i = 0; i < x_sorted.size(); i++) {
        if (i+1 >= x_sorted.size()) {
            // No next vertex to compare with
            break;
        } else if (x_sorted[i].x != x_sorted[i+1].x) {
            // Vertical chords have to have same y coordinate
            continue;
        } else if (x_sorted[i].type != CONCAVE && x_sorted[i+1].type != CONCAVE) {
            // Chords are constructed between two vertices where at least one has to be concave
            continue;
        } else if (x_sorted[i].type == CONVEX || x_sorted[i+1].type == CONVEX) {
            // convex vertices cannot have chords
            continue;
        } else if (AreNeighbors(x_sorted[i], x_sorted[i+1])) {
            // Chords are not constructed between two vertices that are neighbors
            continue;
        } else if (vertex_map_.find(std::make_pair(std::make_pair(x_sorted[i].x, x_sorted[i].y), std::make_pair(x_sorted[i+1].x, x_sorted[i+1].y))) != vertex_map_.end()) {
            continue;
        }

        //std::cout << "Adding vertical chord between vertices (" << x_sorted[i].x << ", " << x_sorted[i].y << ") and (" << x_sorted[i+1].x << ", " << x_sorted[i+1].y << ")\n";
        ver_chords.push_back(std::make_pair(x_sorted[i], x_sorted[i+1]));
        ver_chords_id.push_back(std::make_pair(x_sorted[i].id, x_sorted[i+1].id));
        vertex_map_.insert(std::make_pair(std::make_pair(x_sorted[i].x, x_sorted[i].y), std::make_pair(x_sorted[i+1].x, x_sorted[i+1].y)));
        i += 1; // skip the next vertex as it is already paired
    }


    // Obtain horizontal chords
    std::vector<vertex> y_sorted = vertices_;
    std::sort(y_sorted.begin(), y_sorted.end(), CompVertexOnY);

    for (size_t i = 0; i < y_sorted.size(); i++) {
        if (i+1 >= y_sorted.size()) {
            // No next vertex to compare with
            break;
        } else if (y_sorted[i].y != y_sorted[i+1].y) {
            // Vertical chords have to have same y coordinate
            continue;
        } else if (y_sorted[i].type != CONCAVE && y_sorted[i+1].type != CONCAVE) {
            // Chords are constructed between two vertices where at least one has to be concave
            continue;
        } else if (y_sorted[i].type == CONVEX || y_sorted[i+1].type == CONVEX) {
            continue;
        } else if (AreNeighbors(y_sorted[i], y_sorted[i+1])) {
            // Chords are not constructed between two vertices that are neighbors
            continue;
        } else if (vertex_map_.find(std::make_pair(std::make_pair(y_sorted[i].x, y_sorted[i].y), std::make_pair(y_sorted[i+1].x, y_sorted[i+1].y))) != vertex_map_.end()) {
            continue;
        }
        //std::cout << "Adding horizontal chord between vertices (" << y_sorted[i].x << ", " << y_sorted[i].y << ") and (" << y_sorted[i+1].x << ", " << y_sorted[i+1].y << ")\n";
        hor_chords.push_back(std::make_pair(y_sorted[i], y_sorted[i+1]));
        hor_chords_id.push_back(std::make_pair(y_sorted[i].id, y_sorted[i+1].id));
        vertex_map_.insert(std::make_pair(std::make_pair(y_sorted[i].x, y_sorted[i].y), std::make_pair(y_sorted[i+1].x, y_sorted[i+1].y)));
        i += 1; // skip the next vertex as it is already paired
    }

    // Take minimum
    if (hor_chords.size() < ver_chords.size()) {
        min_chords_ = hor_chords;
        min_chords_id_ = hor_chords_id;
    } else {
        min_chords_ =  ver_chords;
        min_chords_id_ = ver_chords_id;
    }
    return true; 
}


std::vector<vertex> Partitioner::ReduceVertices(std::vector<vertex> vertices) {
    std::vector<vertex> reduced_vertices;
    reduced_vertices.reserve(vertices.size() / 2);

    for (size_t i = 0; i < vertices.size(); i++) {
        if (vertices[i].type != STRAIGHT) {
            reduced_vertices.push_back(vertices[i]);
        }
    }

    return reduced_vertices;
}

std::vector<rect> Partitioner::CornerIdentification(std::vector<vertex> vertices) {
    std::vector<vertex> upper_left, lower_left, lower_right/*,upper_right*/;
    upper_left.reserve(vertices.size() / 4); 
    // upper_right.reserve(vertices.size() / 4);  // upper right not used so commented out for memory/cpu saving
    lower_left.reserve(vertices.size() / 4); 
    lower_right.reserve(vertices.size() / 4);

    for (auto v : vertices) {
        if (v.type == CONVEX) {
            if (v.dir == RIGHT) {
                lower_right.push_back(v);
            } else if (v.dir == DOWN) {
                lower_left.push_back(v);
            } else if (v.dir == LEFT) {
                upper_left.push_back(v);
            } else if (v.dir == UP) {
                // upper_right.push_back(v);
            } else {
                partition_success_ = false;
#ifdef THROW_EXCEPTIONS
                throw std::runtime_error("CornerIdentification fails for vertex with unknown direction");
#endif
#ifdef PRINT_EXCEPTIONS
                std::cout << "CornerIdentification fails for vertex with unknown direction" << std::endl;
#endif
            }
        } 
    }

    for (auto chord : min_chords_id_) {
        vertex v = vertices_[chord.first];
        vertex w = vertices_[chord.second];

        bool vertical;
        if (v.y == w.y) {
            vertical = false;
        } else if (v.x == w.x) {
            vertical = true;
        }

        if (v.dir == RIGHT) {
            if (w.dir == RIGHT) {
                partition_success_ = false;
#ifdef THROW_EXCEPTIONS
                throw std::runtime_error("Chord of vertices with both RIGHT dir");
#endif
#ifdef PRINT_EXCEPTIONS
                std::cout << "Chord of vertices with both RIGHT dir" << std::endl;
#endif
            } else if (w.dir == DOWN && v.type == CONCAVE && w.type == CONCAVE) {
                lower_right.push_back(v);
                //upper_right.push_back(w);
            } else if (w.dir == DOWN && v.type == STRAIGHT && w.type == CONCAVE) {
                lower_left.push_back(v);
                lower_right.push_back(v);
                //upper_right.push_back(w);
            } else if (w.dir == LEFT && v.type == CONCAVE && w.type == CONCAVE && !vertical) {
                upper_left.push_back(v); 
                lower_right.push_back(w);
            } else if (w.dir == LEFT && v.type == CONCAVE && w.type == CONCAVE && vertical) {
                lower_right.push_back(v);
                upper_left.push_back(w);
            } else if (w.dir == LEFT && v.type == CONCAVE && w.type == STRAIGHT) {
                lower_right.push_back(v);
                //upper_right.push_back(w);
                upper_left.push_back(w);
            } else if (w.dir == LEFT && v.type == STRAIGHT && w.type == CONCAVE) {
                lower_left.push_back(v);
                lower_right.push_back(v);
                upper_left.push_back(w);
            } else if (w.dir == UP && v.type == CONCAVE && w.type == CONCAVE) {
                upper_left.push_back(v);
                //upper_right.push_back(w);
            } else if (w.dir == UP && v.type == CONCAVE && w.type == STRAIGHT) {
                upper_left.push_back(v);
                //upper_right.push_back(w);
                lower_right.push_back(w); 
            } else {
                partition_success_ = false;
#ifdef THROW_EXCEPTIONS
                throw std::runtime_error("Unknown case for dir1 " + std::to_string(v.dir) + " type " + std::to_string(v.type) + " and dir2 " + std::to_string(w.dir) + " type " + std::to_string(w.type));
#endif
#ifdef PRINT_EXCEPTIONS
                std::cout << "Unknown case for dir1 " << v.dir << " type " << v.type << " and dir2 " << w.dir << " type " << w.type << std::endl;
#endif
            }
        } else if (v.dir == DOWN) {
            if (w.dir == RIGHT  && v.type == CONCAVE && w.type == CONCAVE) {
                //upper_right.push_back(v);
                lower_right.push_back(w);
            } else if (w.dir == RIGHT  && v.type == CONCAVE && w.type == STRAIGHT) {
                //upper_right.push_back(v);
                lower_right.push_back(w);
                lower_left.push_back(w);
            } else if (w.dir == DOWN) {
                partition_success_ = false;
#ifdef THROW_EXCEPTIONS
                throw std::runtime_error("Chord of vertices with both DOWN dir");
#endif
#ifdef PRINT_EXCEPTIONS
                std::cout << "Chord of vertices with both DOWN dir" << std::endl;
#endif
            } else if (w.dir == LEFT && v.type == CONCAVE && w.type == CONCAVE) {
                lower_left.push_back(v);
                lower_right.push_back(w);
            } else if (w.dir == LEFT && v.type == STRAIGHT && w.type == CONCAVE) {
                lower_left.push_back(v);
                upper_left.push_back(v);
                lower_right.push_back(w);
            } else if (w.dir == UP && v.type == CONCAVE && w.type == CONCAVE) {
                lower_left.push_back(v);
                //upper_right.push_back(w);
            } else if (w.dir == UP && v.type == CONCAVE && w.type == STRAIGHT) {
                lower_left.push_back(v);
                //upper_right.push_back(w);
                lower_right.push_back(w);
            } else if (w.dir == UP && v.type == STRAIGHT && w.type == CONCAVE) {
                upper_left.push_back(v);
                lower_left.push_back(v);
                //upper_right.push_back(w);
            } else {
                partition_success_ = false;
#ifdef THROW_EXCEPTIONS
                throw std::runtime_error("Unknown case for dir1 " + std::to_string(v.dir) + " type " + std::to_string(v.type) + " and dir2 " + std::to_string(w.dir) + " type " + std::to_string(w.type));
#endif
#ifdef PRINT_EXCEPTIONS
                std::cout << "Unknown case for dir1 " << v.dir << " type " << v.type << " and dir2 " << w.dir << " type " << w.type << std::endl;
#endif
            }
        } else if (v.dir == LEFT) {
            if (w.dir == RIGHT && v.type == CONCAVE && w.type == CONCAVE && !vertical) {
                lower_right.push_back(v);
                upper_left.push_back(w);
            } else if (w.dir == RIGHT && v.type == CONCAVE && w.type == CONCAVE && vertical) {
                upper_left.push_back(v);
                lower_right.push_back(w);
            } else if (w.dir == RIGHT && v.type == CONCAVE && w.type == STRAIGHT) {
                lower_right.push_back(v);
                lower_left.push_back(v);
                //upper_right.push_back(w);
            } else if (w.dir == RIGHT && v.type == STRAIGHT && w.type == CONCAVE) {
                //upper_right.push_back(v);
                upper_left.push_back(v);
                lower_right.push_back(w);
            } else if (w.dir == DOWN && v.type == CONCAVE && w.type == STRAIGHT) {
                lower_right.push_back(v);
                lower_left.push_back(w);
                upper_left.push_back(w);
            } else if (w.dir == LEFT) {
                partition_success_ = false;
#ifdef THROW_EXCEPTIONS
                throw std::runtime_error("Chord of vertices with both LEFT dir");
#endif
#ifdef PRINT_EXCEPTIONS
                std::cout << "Chord of vertices with both LEFT dir" << std::endl;
#endif
            } else if (w.dir == UP && v.type == CONCAVE && w.type == CONCAVE) {
                upper_left.push_back(v);
                lower_left.push_back(w);
            } else if (w.dir == UP && v.type == STRAIGHT && w.type == CONCAVE) {
                //upper_right.push_back(v);
                upper_left.push_back(v);
                lower_right.push_back(w);
            } else {
                partition_success_ = false;
#ifdef THROW_EXCEPTIONS
                throw std::runtime_error("Unknown case for dir1 " + std::to_string(v.dir) + " type " + std::to_string(v.type) + " and dir2 " + std::to_string(w.dir) + " type " + std::to_string(w.type));
#endif
#ifdef PRINT_EXCEPTIONS
                std::cout << "Unknown case for dir1 " << v.dir << " type " << v.type << " and dir2 " << w.dir << " type " << w.type << std::endl;
#endif
            }
        } else if (v.dir == UP) {
            if (w.dir == RIGHT && v.type == STRAIGHT && w.type == CONCAVE) {
                //upper_right.push_back(v);
                lower_right.push_back(v);
                upper_left.push_back(w);
            } else if (w.dir == DOWN && v.type == CONCAVE && w.type == CONCAVE) {
                lower_left.push_back(v);
                //upper_right.push_back(w);
            } else if (w.dir == DOWN && v.type == CONCAVE && w.type == STRAIGHT) {
                //upper_right.push_back(v);
                lower_left.push_back(w);
                upper_left.push_back(w);
            } else if (w.dir == DOWN && v.type == STRAIGHT && w.type == CONCAVE) {
                lower_right.push_back(v);
                //upper_right.push_back(v);
                lower_left.push_back(w);
            } else if (w.dir == LEFT && v.type == CONCAVE && w.type == CONCAVE) {
                lower_left.push_back(v);
                upper_left.push_back(w);
            } else if (w.dir == LEFT && v.type == CONCAVE && w.type == STRAIGHT) {
                lower_left.push_back(v);
                upper_left.push_back(w);
                //upper_right.push_back(w);
            } else if (w.dir == UP) {
                partition_success_ = false;
#ifdef THROW_EXCEPTIONS
                throw std::runtime_error("Chord of vertices with both UP dir");
#endif
#ifdef PRINT_EXCEPTIONS
                std::cout << "Chord of vertices with both UP dir" << std::endl;
#endif
            } else {
                partition_success_ = false;
#ifdef THROW_EXCEPTIONS
                throw std::runtime_error("Unknown case for dir1 " + std::to_string(v.dir) + " type " + std::to_string(v.type) + " and dir2 " + std::to_string(w.dir) + " type " + std::to_string(w.type));
#endif
#ifdef PRINT_EXCEPTIONS
                std::cout << "Unknown case for dir1 " << v.dir << " type " << v.type << " and dir2 " << w.dir << " type " << w.type << std::endl;
#endif
            }
        }
    }

    // Sort upper left points on x coordinate, then on y
    // Sort lower right by y-coordinate, then on x
    std::sort(upper_left.begin(), upper_left.end(), CompVertexOnX);
    std::sort(lower_right.begin(), lower_right.end(), CompVertexOnY);

    // Rectangle construction for point p in lower left 
    // for p in lower left:
    //      Find the point in upper_left closest to p but greater than it, priotizing x-coord over y-coordinate, denote u1
    //      Find the point in lower_left closest to p but greater than it, prioritizing y-coord over x, denote u2
    //      Add rectangle with upper-left u1, lower-right u2 point to R
    std::vector<rect> rectangles;
    rectangles.reserve(lower_left.size());
    for (auto v : lower_left) {
        auto ul = std::upper_bound(upper_left.begin(), upper_left.end(), v, CompVertexOnX);
        auto lr = std::upper_bound(lower_right.begin(), lower_right.end(), v, CompVertexOnY);

        if (ul == upper_left.end() || lr == lower_right.end()) {
            std::cout << "Cannot find successive point for vertex " << v.id << ", x: " << v.x << ", y: " << v.y << std::endl;
            partition_success_ = false;
#ifdef THROW_EXCEPTIONS
            throw std::runtime_error("Rectangle extraction failed");
#endif
            continue;
        }

        // Obtain lower left corner and upper right from upper left and lower right
        rect r = {v.x, v.y, (*lr).x, (*ul).y}; //{(*ul).x, (*lr).y, (*lr).x, (*ul).y};
        rectangles.push_back(r);
    }

    if (lower_left.size() != lower_right.size() || lower_left.size() != upper_left.size()) {
#ifdef PRINT_EXCEPTIONS
        std::cout << "Partition fail on non-equal corners" << std::endl;
#endif
        partition_success_ = false;
    } else if (lower_left.size() != rectangles.size()) {
#ifdef PRINT_EXCEPTIONS
        std::cout << "Partition fail on more corners than rectangles" << std::endl;
#endif
        partition_success_ = false;
    }

    return rectangles;
}


void Partitioner::PrintRectangles(std::vector<rect> rectangles) {
    std::cout << rectangles.size() << " rectangles:\n";
    for (rect r : rectangles) {
        std::cout << "(" << r.x1 << ", " << r.y1 << "; " << r.x2 << ", " << r.y2 << ") - ";
    }
    std::cout << std::endl;
}