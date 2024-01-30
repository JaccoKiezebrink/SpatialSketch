#include "RectilinearPartitioner.h"

#include <iostream>
#include <math.h>
#include <algorithm>
#include <list>
#include <string.h>
#include <stdexcept>
#include <unistd.h>

/* OLD RECTILINEAR POLYGON PARTITIONING CODE, CORRECT AND MORE EFFICIENT CODE IN PARTITIONER.CPP */


//#define DEBUG  // Uncomment if print statements are desired

/* TODO:
    * If a horizontal chord is picked and there needs to be a vertical chord added to ensure rectangles, this new chord intersects with the horizontal chord,
      then there is no mechanism that lets the chord be drawn only up until the horizontal chord.
*/

RectilinearPartitioner::RectilinearPartitioner(std::vector<std::pair<float, float>> vertices) {
    vertices_ = vertices;
    vertex_to_chord_.resize(vertices_.size());
    vertex_to_matched_chord_.resize(vertices_.size());

    if (!ComputeVertexType()) {
        std::cout << "Vertex type computation went wrong\n" << std::endl;
        //throw std::exception();
        success_ = false;
    }
}

RectilinearPartitioner::~RectilinearPartitioner() {
    vertices_.clear(); 
    convex_vertices_id_.clear();
    concave_vertices_id_.clear();
    straight_vertices_id_.clear();
    concave_vertices_.clear();
    hor_chords_.clear();
    ver_chords_.clear();
    matching_hor_chords_.clear();
    matching_ver_chords_.clear();
    vertex_to_chord_.clear();
    vertex_to_matched_chord_.clear();
    independent_chords_.clear();
    final_chord_set_.clear();
}

std::vector<rect> RectilinearPartitioner::ReturnRect() {
        std::vector<std::pair<float, float>> corners;
        std::vector<rect> result;
        for (int i = 0; i < (int) convex_vertices_id_.size(); i++) {
            corners.push_back(vertices_[convex_vertices_id_[i]]);
        }
        result.push_back(CornersToRectangle(corners));

#ifdef DEBUG
        std::cout << "Shape only has four convex corner, thus rectangle is assumed." << std::endl;
#endif

        return result;
}

std::vector<rect> RectilinearPartitioner::ComputePartitioning(int method) {
    if (method == 1) {
        return ComputePartitioningMinChord();
    } else {
        return ComputePartitioningMNC();
    }
}

std::vector<rect> RectilinearPartitioner::ComputePartitioningMNC() {
    // If range is rectangular
    if (convex_vertices_id_.size() == 4 && concave_vertices_id_.size() == 0) {
        return ReturnRect();
    }
    ComputeChords();
    MaximumIndependentSet();
    FinalizeChordSet();
    return GetRectangleCoordinates();
}


std::vector<rect> RectilinearPartitioner::ComputePartitioningMinChord() {
    // If range is rectangular
    if (convex_vertices_id_.size() == 4 && concave_vertices_id_.size() == 0) {
        return ReturnRect();
    }
    ComputeChords();

    /*std::cout << "Vertical chords: " << std::endl;
    for (auto chord : ver_chords_) {
        std::cout << vertices_[chord.first].first << "," << vertices_[chord.first].second << vertices_[chord.second].first << "," << vertices_[chord.second].second << std::endl;
    }
    std::cout << "Horizontal chords: " << std::endl;
    for (auto chord : hor_chords_) {
        std::cout << vertices_[chord.first].first << "," << vertices_[chord.first].second << vertices_[chord.second].first << "," << vertices_[chord.second].second << std::endl;
    }*/

    /*if (hor_chords_.size() <= ver_chords_.size()) {
        for (int i = 0; i < (int) hor_chords_.size(); i++) {
            final_chord_set_.push_back(i);
        }
    } else {
        for (int i = 0; i < (int) ver_chords_.size(); i++) {
            final_chord_set_.push_back(i + hor_chords_.size());
        }
    }*/

    FinalizeChordSet();
    return GetRectangleCoordinates();
}

bool RectilinearPartitioner::ComputeVertexType() {
// Given vector of vertices (x, y) sorted in anti-clockwise order, compute for each vertex if it is a concave vertex, that is, if it creates an internal angle of 270 degrees.
    int length = vertices_.size();
    for (int i = 1; i < (int) vertices_.size() + 1; i++) {
        std::pair<float, float> v1, v2, v3, diff1, diff2; 
        Direction dir1, dir2;  // direction of v1-v2 and v2-v3
        v1 = vertices_[i%length];
        v2 = vertices_[(i+1)%length];
        v3 = vertices_[(i+2)%length];

        //std::cout << "Coordinates: (" << v2.first << ", " << v2.second << ") ";

        diff1.first = v2.first - v1.first;
        diff1.second = v2.second - v1.second;
        diff2.first = v3.first - v2.first;
        diff2.second = v3.second - v2.second;

        if (diff1.first == 0 && diff1.second == 1) {
            dir1 = UP;
        } else if (diff1.first == 0 && diff1.second == -1) {
            dir1 = DOWN;
        } else if (diff1.first == 1 && diff1.second == 0) {
            dir1 = RIGHT;
        } else if (diff1.first == -1 && diff1.second == 0) {
            dir1 = LEFT;
        } else {
            std::cout << "!Error: (" << diff1.first << ", " << diff1.second << ") is not a valid direction" << std::endl;
            std::cout << "Coordinates: (" << v1.first << ", " << v1.second << "), (" << v2.first << ", " << v2.second << "), (" << v3.first << ", " << v3.second << ")" << std::endl;
            return false;
        }

        if (diff2.first == 0 && diff2.second == 1) {
            dir2 = UP;
        } else if (diff2.first == 0 && diff2.second == -1) {
            dir2 = DOWN;
        } else if (diff2.first == 1 && diff2.second == 0) {
            dir2 = RIGHT;
        } else if (diff2.first == -1 && diff2.second == 0) {
            dir2 = LEFT;
        } else {
            std::cout << "!Error: (" << diff2.first << ", " << diff2.second << ") is not a valid direction" << std::endl;
            std::cout << "Coordinates: (" << v1.first << ", " << v1.second << "), (" << v2.first << ", " << v2.second << "), (" << v3.first << ", " << v3.second << ")" << std::endl;
            return false;
        }

        if (dir1 == dir2) {
            straight_vertices_id_.push_back((i+1)%length);
#ifdef DEBUG
            std::cout << (i+1)%length << " straight" << std::endl;
#endif
        } else if ((dir1 == RIGHT && dir2 == DOWN) || (dir1 == LEFT && dir2 == UP) || (dir1 == UP && dir2 == RIGHT) || (dir1 == DOWN && dir2 == LEFT)) {
            concave_vertices_id_.push_back((i+1)%length);
            concave_vertices_.push_back(std::make_tuple(v2.first, v2.second, (i+1)%length));
#ifdef DEBUG
            std::cout << (i+1)%length << " concave" << std::endl;
#endif
        } else if ((dir1 == RIGHT && dir2 == UP) || (dir1 == LEFT && dir2 == DOWN) || (dir1 == UP && dir2 == LEFT) || (dir1 == DOWN && dir2 == RIGHT)) {
            convex_vertices_id_.push_back((i+1)%length);
#ifdef DEBUG
            std::cout << (i+1)%length << " convex" << std::endl;
#endif
        } else {
            // Optional assume straight: 
            straight_vertices_id_.push_back((i+1)%length);
            //concave_vertices_id_.push_back((i+1)%length);
            //concave_vertices_.push_back(std::make_tuple(v2.first, v2.second, (i+1)%length));
            //std::cout << "!Error: (" << dir1 << ", " << dir2 << ") by vertices " << (i)%length << ", " << (i+1)%length << ", " << (i+2)%length << " does not create valid angle" << std::endl;
            //std::cout << "Coordinates: (" << v1.first << ", " << v1.second << "), (" << v2.first << ", " << v2.second << "), (" << v3.first << ", " << v3.second << ")" << std::endl;
            //return false;
        }
    }
    return true;
}


void RectilinearPartitioner::ComputeChords() {
    if (hor_chords_.size() > 0 || ver_chords_.size() > 0) {
        hor_chords_.clear();
        ver_chords_.clear();
        matching_hor_chords_.clear();
        matching_ver_chords_.clear();
    }

    // vertex to vertical chord mapping, later to be merge in combined mapping
    std::vector<int> vertex_to_ver_chord(vertices_.size(), -1);

    // Store vertices sorted on both x as y coordinates locally
    std::vector<std::pair<float, float>> vertices_sorted_x = vertices_;
    std::vector<std::pair<float, float>> vertices_sorted_y = vertices_;
    std::sort(vertices_sorted_x.begin(), vertices_sorted_x.end(), CompByFirstFloat);
    std::sort(vertices_sorted_y.begin(), vertices_sorted_y.end(), CompBySecondFloat);

    // sort concave vertices on x
    std::vector<std::tuple<float, float, int>> concave_vertices_sorted_x = concave_vertices_;
    std::vector<std::tuple<float, float, int>> concave_vertices_sorted_y = concave_vertices_;
    std::sort(concave_vertices_sorted_x.begin(), concave_vertices_sorted_x.end(), CompByFirstTupleFloat);
    std::sort(concave_vertices_sorted_y.begin(), concave_vertices_sorted_y.end(), CompBySecondTupleFloat);

    // extract chords between two concave vertices
    for (int i = 0; i < (int) concave_vertices_id_.size(); i++) {
        std::pair<float, float> v1, v2;
        v1 = vertices_[concave_vertices_id_[i]];
        
        bool intersect = false;

        // Binary search other concave vertices with the same x value
        auto vertices_x = std::equal_range(concave_vertices_sorted_x.begin(), concave_vertices_sorted_x.end(), v1.first, CompFirstFloatTuple{});
        for (auto it = vertices_x.first; it != vertices_x.second; it++) {
            int j = std::get<2>(*it);
            // Ignore identical vertices or neighboring vertices
            if (j <= concave_vertices_id_[i] || std::abs(j - concave_vertices_id_[i]) == 1 || std::abs(j - concave_vertices_id_[i]) == vertices_.size() - 1) {
                continue;
            }

            v2 = std::make_pair(std::get<0>(*it), std::get<1>(*it));

            // if x coordinates are identical
            intersect = false;
            auto range = std::equal_range(vertices_sorted_x.begin(), vertices_sorted_x.end(), v1.first, CompFirstFloatPair{});
            for (auto vec = range.first; vec != range.second; vec++) {
                // check if vertex has same x coordinate as well and lies between v1 and v2 on the y-axis
                if (vec->first == v1.first && ((v1.second > vec->second && v2.second < vec->second) || (v1.second < vec->second && v2.second > vec->second))) {
                    intersect = true;
                    break;
                }
            }
            if (!intersect) { // middles.size() <= 0
                vertex_to_ver_chord[concave_vertices_id_[i]] = ver_chords_.size();
                vertex_to_ver_chord[j] = ver_chords_.size();

                int min = std::min(concave_vertices_id_[i], j);
                int max = std::max(concave_vertices_id_[i], j);
                ver_chords_.push_back(std::make_pair(min, max));
                matching_ver_chords_.push_back(std::make_pair(min, max));
                break;
            }
        }
    }

    for (int i = 0; i < (int) concave_vertices_id_.size(); i++) {
        std::pair<float, float> v1, v2;
        v1 = vertices_[concave_vertices_id_[i]];

        bool intersect = false;

        auto vertices_y = std::equal_range(concave_vertices_sorted_y.begin(), concave_vertices_sorted_y.end(), v1.second, CompSecondFloatTuple{});
        for (auto it = vertices_y.first; it != vertices_y.second; it++) {
            int j = std::get<2>(*it);
            // Ignore identical vertices and neighboring vertices
            if (j <= concave_vertices_id_[i] || std::abs(j - concave_vertices_id_[i]) == 1 || std::abs(j - concave_vertices_id_[i]) == vertices_.size() - 1) {
                continue;
            }

            v2 = std::make_pair(std::get<0>(*it), std::get<1>(*it));

            // Do the same but on y-axis
            intersect = false;
            auto range = std::equal_range(vertices_sorted_y.begin(), vertices_sorted_y.end(), v1.second, CompSecondFloatPair{});
            for (auto vec = range.first; vec != range.second; vec++) {
                // check if vertex has same x coordinate as well and lies between v1 and v2 on the y-axis
                if (vec->second == v1.second && ((v1.first > vec->first && v2.first < vec->first) || (v1.first < vec->first && v2.first > vec->first))) {
                    intersect = true;
                    break;
                }
            }
            if (!intersect) { // middles.size() <= 0
                vertex_to_chord_[concave_vertices_id_[i]].push_back(hor_chords_.size());
                vertex_to_chord_[j].push_back(hor_chords_.size());

                int min = std::min(concave_vertices_id_[i], j);
                int max = std::max(concave_vertices_id_[i], j);
                hor_chords_.push_back(std::make_pair(min, max));
                matching_hor_chords_.push_back(std::make_pair(min, max));
                break;
            }
        }
    }

    // Setup vertex to matched chords mapping
    for (int i = 0; i < (int) matching_hor_chords_.size(); i++) {
        std::pair<int, int> chord = matching_hor_chords_[i];
        vertex_to_matched_chord_[chord.first].push_back(i);
        vertex_to_matched_chord_[chord.second].push_back(i);
    }

    for (int i = i; i < (int) matching_ver_chords_.size(); i++) {
        vertex_to_matched_chord_[matching_ver_chords_[i].first].push_back(i + matching_hor_chords_.size());
        vertex_to_matched_chord_[matching_ver_chords_[i].second].push_back(i + matching_hor_chords_.size());
    }


    // Extract chords between straight_vertices_and concave vertices
    for (int i = 0; i < (int) straight_vertices_id_.size(); i++) {
        std::pair<float, float> v1, v2;
        v1 = vertices_[straight_vertices_id_[i]];
        bool intersect = false;

        // If the straight vertex lays vertically, we only have to find other vertices with the same x coordinate
        Orientation orientation = DetermineOrientation(straight_vertices_id_[i], (straight_vertices_id_[i]+1) % vertices_.size());
        if (orientation == VERTICAL) {
            // Find concave vertices with same x coordinate
            auto vertices_x = std::equal_range(concave_vertices_sorted_x.begin(), concave_vertices_sorted_x.end(), v1.first, CompFirstFloatTuple{});
            for (auto it = vertices_x.first; it != vertices_x.second; it++) {
                int j = std::get<2>(*it);
                v2 = std::make_pair(std::get<0>(*it), std::get<1>(*it));

                // Skip neighboring vertices
                if (std::abs(concave_vertices_id_[j] - straight_vertices_id_[i]) == 1 || std::abs(concave_vertices_id_[j] - straight_vertices_id_[i]) == vertices_.size() - 1) {
                    continue;
                }

                // Find all other vertices and compare whether they lay within v1 and v1, if so set intersect to true
                auto range = std::equal_range(vertices_sorted_x.begin(), vertices_sorted_x.end(), v1.first, CompFirstFloatPair{});
                intersect = false;
                for (auto vec = range.first; vec != range.second; vec++) {
                    // check if vertex has same x coordinate as well and lies between v1 and v2 on the y-axis
                    if (vec->first == v1.first && ((v1.second > vec->second && v2.second < vec->second) || (v1.second < vec->second && v2.second > vec->second))) {
                        intersect = true;
                        break;
                    }
                }

                // If no intersecting vertex was found, add a chord between v1 and v2
                if (!intersect) {
                    vertex_to_ver_chord[straight_vertices_id_[i]] = ver_chords_.size();
                    vertex_to_ver_chord[j] = ver_chords_.size();

                    int min = std::min(straight_vertices_id_[i], j);
                    int max = std::max(straight_vertices_id_[i], j);
                    ver_chords_.push_back(std::make_pair(min, max));
                    break;
                }
            }
        } else if (orientation == HORIZONTAL) {  // straight vertex lays horizontally thus we only have to find other vertices with the same y coordinate 
            // Find vertices with equivalent y coordinat
            auto vertices_y = std::equal_range(concave_vertices_sorted_y.begin(), concave_vertices_sorted_y.end(), v1.second, CompSecondFloatTuple{});
            for (auto it = vertices_y.first; it != vertices_y.second; it++) {
                int j = std::get<2>(*it);
                v2 = std::make_pair(std::get<0>(*it), std::get<1>(*it));

                // Check if vertex are neighbors
                if (std::abs(j - straight_vertices_id_[i]) == 1 || std::abs(j - straight_vertices_id_[i]) == vertices_.size() - 1) {
                    continue;
                }

                // obtain all vertices on this y coordinate
                auto range = std::equal_range(vertices_sorted_y.begin(), vertices_sorted_y.end(), v1.second, CompSecondFloatPair{});
                intersect = false;
                for (auto vec = range.first; vec != range.second; vec++) {
                    // check if vertex has same x coordinate as well and lies between v1 and v2 on the y-axis
                    if (vec->second == v1.second && ((v1.first > vec->first && v2.first < vec->first) || (v1.first < vec->first && v2.first > vec->first))) {
                        intersect = true;
                    }
                }
                // If no vertex is between v1 and v2
                if (!intersect) {
                    vertex_to_chord_[straight_vertices_id_[i]].push_back(hor_chords_.size());
                    vertex_to_chord_[j].push_back(hor_chords_.size());

                    int min = std::min(straight_vertices_id_[i], j);
                    int max = std::max(straight_vertices_id_[i], j);
                    hor_chords_.push_back(std::make_pair(min, max));
                    break;
                }
            }
        } else {
            std::pair<float, float> first = vertices_[straight_vertices_id_[i]];
            std::pair<float, float> second = vertices_[(straight_vertices_id_[i] + 1) % vertices_.size()];
            std::cout << "Orientation of vertices " << straight_vertices_id_[i] << " (" << first.first << ", " << first.second << " and " << (straight_vertices_id_[i]+1) % vertices_.size() << " (" << second.first << ", " << second.second << ") could not be determined" << std::endl;
        }
    }

    // Append vertex to vertical chord mapping in combined mapping
    int H = hor_chords_.size();
    for (int i = 0; i < (int) vertex_to_ver_chord.size(); i++) {
        if (vertex_to_ver_chord[i] != -1) {
            vertex_to_chord_[i].push_back((vertex_to_ver_chord[i] + H));
        }
    }

#ifdef DEBUG
    std::cout << "\nhor chords" << std::endl;
    for (int i = 0; i < (int) hor_chords_.size(); i++) {
        //std::cout << hor_chords_[i].first << ", " << hor_chords_[i].second << std::endl;
        std::cout << i << ": (" << vertices_[hor_chords_[i].first].first << ", " << vertices_[hor_chords_[i].first].second << ") - (" << vertices_[hor_chords_[i].second].first << ", " << vertices_[hor_chords_[i].second].second << ")";
        for (int j = 0; j < (int) matching_hor_chords_.size(); j++) {
            if (matching_hor_chords_[j] == hor_chords_[i]) {
                std::cout << " (in matching)";
            }
        }
        std::cout << std::endl;
    }
    std::cout << "ver chords" << std::endl;
    for (int i = 0; i < (int) ver_chords_.size(); i++) {
        //std::cout << ver_chords_[i].first << ", " << ver_chords_[i].second << std::endl;
        std::cout << i + hor_chords_.size() << ": (" << vertices_[ver_chords_[i].first].first << ", " << vertices_[ver_chords_[i].first].second << ") - (" << vertices_[ver_chords_[i].second].first << ", " << vertices_[ver_chords_[i].second].second << ")";
        for (int j = 0; j < (int) matching_ver_chords_.size(); j++) {
            if (matching_ver_chords_[j] == ver_chords_[i]) {
                std::cout << " (in matching)";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
#endif
}


void RectilinearPartitioner::MaximumIndependentSet() {
    int H = matching_hor_chords_.size();
    int V = matching_ver_chords_.size();
    bool **graph = new bool*[H];
    for (int i = 0; i < H; i++) {
        graph[i] = new bool[V];
    }

    // For all chords, verify if they intersect and if so add them to the graph
    for (int i = 0; i < (int) matching_hor_chords_.size(); i++) {
        for (int j = 0; j < (int) matching_ver_chords_.size(); j++) {
            std::pair<float, float> hor_v1, hor_v2, ver_v1, ver_v2;
            hor_v1 = vertices_[matching_hor_chords_[i].first];
            hor_v2 = vertices_[matching_hor_chords_[i].second];
            ver_v1 = vertices_[matching_ver_chords_[j].first];
            ver_v2 = vertices_[matching_ver_chords_[j].second];

            // Check if x of vertical chord lays in range of horizontal chord and y of horizontal chord lays in range of vertical chord
            if ((ver_v1.first >= std::min(hor_v1.first, hor_v2.first) && ver_v1.first <= std::max(hor_v1.first, hor_v2.first))
                && (hor_v1.second >= std::min(ver_v1.second, ver_v2.second) && hor_v1.second <= std::max(ver_v1.second, ver_v2.second))) {
                graph[i][j] = true;
            } else {
                graph[i][j] = false;
            }
        }
    }

#ifdef DEBUG
    // Print graph
    std::cout << "Matching graph" << std::endl;
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < V; j++) {
            std::cout << graph[i][j] << " ";
        }
        std::cout << std::endl;
    }
#endif

    // Following bipartite matching code retrieved from: https://www.geeksforgeeks.org/maximum-bipartite-matching/
    // array that keeps track of vertical chords matched, with -1 indicating no match and some value indicating the horizontal chord
    int matchH[H], matchV[V];
    memset(matchH, -1, sizeof(matchH));
    memset(matchV, -1, sizeof(matchV));

    int nr_of_matches = 0;
    for (int h = 0; h < H; h++) {
        // Mark all nodes as not seen for next chord.
        bool seen[V];
        memset(seen, 0, sizeof(seen));
 
        // Find if the applicant 'u' can get a job
        if (BipartiteMatch(graph, h, seen, matchH, matchV)) {
            nr_of_matches++;
        }
    }

    for (int i = 0; i < V; i++) {
        //std::cout << "matchV[" << i << "] (" << matchV[i] << ") != -1 == " << (matchV[i] != -1) << std::endl;
        if (matchV[i] != -1) {
            //std::cout << "matchH[" << matchV[i] << "] (" << matchH[matchV[i]] << ") = " << i << std::endl;
            matchH[matchV[i]] = i;
        }
        //std::cout << std::endl;
        //std::cout << "Set matchH[" << matchV[i] << "] = " << i << std::endl;
    }

#ifdef DEBUG
    std::cout << nr_of_matches << std::endl;
    for (int i = 0; i < V; i++) {
        if (matchV[i] != -1) {
            std::cout << "v chord " << i + H << " matched with h chord " << matchV[i] << " --> check: " << matchH[matchV[i]] << std::endl;
        } else {
            std::cout << "v chord " << i + H << " not matched" << std::endl;
        }
    }
#endif

    std::vector<int> independent_vertices, unusable_chords;
    std::vector<int> free_vertices;
    for (int i = 0; i < H; i++) {
        // if vertex is unmatched and has no edges
        if (matchH[i] == -1) {
            // then it's definitely a free vertex
            free_vertices.push_back(i);
        }
    }
    
    for (int i = 0; i < V; i++) {
        if (matchV[i] == -1) {
            free_vertices.push_back(i + H);
        }
    }

#ifdef DEBUG
    std::cout << "initially free vertices (" << free_vertices.size() << "): ";
    for (auto v : free_vertices) {
        std::cout << v << ", ";
    }
    std::cout << std::endl;
#endif 

    while (free_vertices.size() != 0 || nr_of_matches) {
        int vertex_u = -1;
        int vertex_v;
        // If there are free vertices, process these
        if (free_vertices.size() > 0) {
            vertex_u = free_vertices.back();
            free_vertices.pop_back();
            independent_vertices.push_back(vertex_u);
#ifdef DEBUG
            std::cout << "Vertex u " << vertex_u << " is free" << std::endl;
#endif
        } else { // otherwise, process a matched vertex
            for (int i = 0; i < H; i++) {
                if (matchH[i] != -1) {
                    vertex_u = i;
                    vertex_v = matchH[i];

                    matchH[i] = -1;
                    matchV[vertex_v] = -1;
                    
                    // If edge has already been removed, continue
                    if (!graph[vertex_u][vertex_v]) {
                        continue;
                    }

                    nr_of_matches--;
                    break;
                }
            }

            if (vertex_u == -1) {
                std::cout << "\n!No vertex_u found, but nr of matches still " << nr_of_matches << std::endl;
                break;
            }

            graph[vertex_u][vertex_v] = false;
            independent_vertices.push_back(vertex_u);
#ifdef DEBUG
            std::cout << "Vertex u " << vertex_u << " picked from matching (vertex v " << vertex_v << ")" << std::endl;
#endif
        }

        // Process vertex u
        if (vertex_u < H) {
            for (int j = 0; j < V; j++) {
                if (graph[vertex_u][j]) {
                    graph[vertex_u][j] = false;
                    if (matchV[j] != -1) {
#ifdef DEBUG
                       std::cout << "Vertex v " << j << " is matched, so chord " << j+H << " becomes unusable and " << matchV[j] << " becomes free" << std::endl;
#endif
                        if (matchV[j] != -1) {
                            free_vertices.push_back(matchV[j]);
                            matchV[j] = -1;
                        }
                        unusable_chords.push_back(j + hor_chords_.size());  // unusable chords stores chords as their actual number not only mathced number
                        matchH[matchV[j]] = -1;
                        nr_of_matches--;
                    }
                }
            }
        } else {
            for (int j = 0; j < H; j++) {
                if (graph[j][vertex_u]) {
                    graph[j][vertex_u] = false;
                    if (matchH[j] != -1) {
#ifdef DEBUG
                        std::cout << "Vertex h " << j << " is matched so it becomes unusable, and " << matchH[j] + H << " becomes free" << std::endl;
#endif
                        free_vertices.push_back(matchH[j] + H);
                        unusable_chords.push_back(j);
                        matchV[matchH[j]] = -1;
                        matchH[j] = -1;
                        nr_of_matches--;
                    }
                }
            }
        }
    }

#ifdef DEBUG
    std::cout << "Independent set " << independent_vertices.size() << std::endl;
    for (int i = 0; i < (int) independent_vertices.size(); i++) {
        std::cout << independent_vertices[i] << ", ";
    }
    std::cout << std::endl << std::endl;
#endif     

   independent_chords_ = independent_vertices;
   std::sort(unusable_chords.begin(), unusable_chords.end());
   unusable_chords_ = unusable_chords;
}

// https://www.geeksforgeeks.org/maximum-bipartite-matching/
// A DFS based recursive function that returns true if a matching for vertex u is possible
bool RectilinearPartitioner::BipartiteMatch(bool **graph, int h, bool seen[], int matchH[], int matchV[]) {
    // Try every chord one by one
    for (int v = 0; v < (int) ver_chords_.size(); v++) {
        // If edge is unvisited
        if (graph[h][v] && !seen[v]) {
            // Mark v as visited
            seen[v] = true;

            // If chord v is not assigned an edge yet or previously assigned has alternate edge available
            // Since v is marked as visited in the above line, matchR[v] in the following recursive call will not get job 'v' again
            if (matchV[v] < 0 || BipartiteMatch(graph, matchV[v], seen, matchH, matchV)) {
                matchV[v] = h;
                //matchH[h] = v;
                //std::cout << "matchV[" << v <<"] = " << h << ", matchH[" << h << "] = " << v << std::endl;
                return true;
            }
        }
    }
    return false;
}


bool RectilinearPartitioner::CheckIntersect(std::pair<int, int> chord1, std::pair<int, int> chord2) {
    Orientation o1 = DetermineOrientation(chord1.first, chord1.second);
    Orientation o2 = DetermineOrientation(chord2.first, chord2.second);

    // If orientation is the same they cannot overlap
    if (o1 == o2) {
        return false;
    }

    std::pair<int, int> c1, c2;
    if (o1 == VERTICAL) {
        c1 = chord1;
        c2 = chord2;   
    } else {
        c1 = chord2;
        c2 = chord1;
    }

    // if x of chord 1 lies in between the two x-values of chord 2
    if ((vertices_[c1.first].first < vertices_[c2.first].first && vertices_[c1.second].first > vertices_[c2.first].first) || 
        (vertices_[c1.first].first > vertices_[c2.first].first && vertices_[c1.second].first < vertices_[c2.first].first)) {
        // then y of chord 2 has to lay in between y values of chord 1
        if ((vertices_[c2.first].second < vertices_[c1.first].second && vertices_[c2.second].second > vertices_[c1.first].second) || 
            (vertices_[c2.first].second > vertices_[c1.first].second && vertices_[c2.second].second < vertices_[c1.first].second)) {
            return true;
        }
    }
    return false;    
}


void RectilinearPartitioner::FinalizeChordSet() {
    std::vector<int> final_chord_set;
    bool has_chord[vertices_.size()] = {false};
    int H = matching_hor_chords_.size();

    // Mark concave vertices with selected chord as visited
    for (int i = 0; i < (int) independent_chords_.size(); i++) {
        if (independent_chords_[i] < H) {
            auto chord = matching_hor_chords_[independent_chords_[i]];
            has_chord[chord.first] = true;
            has_chord[chord.second] = true;
#ifdef DEBUG
            std::cout << "Adding horizontal chord (" << vertices_[chord.first].first << ", " << vertices_[chord.first].second << "), (" << vertices_[chord.second].first << ", " << vertices_[chord.second].second << ") to final set due to matching" << std::endl;
#endif
        } else {
            auto chord = matching_ver_chords_[independent_chords_[i] - H];
            has_chord[chord.first] = true;
            has_chord[chord.second] = true;
#ifdef DEBUG
            std::cout << "Adding vertical chord (" << vertices_[chord.first].first << ", " << vertices_[chord.first].second << "), (" << vertices_[chord.second].first << ", " << vertices_[chord.second].second << ") to final set due to matching" << std::endl;
#endif
        }
    }

    for (int i = 0; i < (int) concave_vertices_id_.size(); i++) {
        // If concave vertex i does not have a chord, find a horizontal chord and add it to the set
        if (!has_chord[concave_vertices_id_[i]]) {
            std::vector<int> options = vertex_to_chord_[concave_vertices_id_[i]];
            std::pair<int, int> final_chord;
            for (int chord : options) {
                // check if chord has been discarded during matching
                if (std::binary_search(unusable_chords_.begin(), unusable_chords_.end(), chord)) {
#ifdef DEBUG
                    std::cout << "ignoring chord " << chord << " because it is unusable" << std::endl;
#endif
                    continue;
                }

                // Check if it already intersects with any chord in the independent_chords_
                bool no_intersect = true;
                if (chord < (int) hor_chords_.size()) {
                    for (int j = 0; j < (int) independent_chords_.size(); j++) {
                        if (independent_chords_[j] >= H) {  // horizontal chords can only intersect with vertical chords
                            if (CheckIntersect(hor_chords_[chord], matching_ver_chords_[independent_chords_[j]])) {
                                no_intersect = false;
#ifdef DEBUG
                                std::cout << "ignoring hchord " << chord << " because it intersects with " << independent_chords_[j] << std::endl;
#endif
                                break;
                            }
                        }
                    }
                } else {
                    for (int j = 0; j < (int) independent_chords_.size(); j++) {
                        if (independent_chords_[j] < H) {  // vertical chords can only intersect with horizontal chords
                            if (CheckIntersect(ver_chords_[chord - hor_chords_.size()], matching_hor_chords_[independent_chords_[j]])) {
                                no_intersect = false;
#ifdef DEBUG
                                std::cout << "ignoring vchord " << chord << " because it intersects with " << independent_chords_[j] << std::endl;
#endif
                                break;
                            }
                        }
                    }
                }

                if (no_intersect) {
                    final_chord_set.push_back(chord);

                    // mark vertices now have a chord
                    if (chord < (int) hor_chords_.size()) {
                        has_chord[hor_chords_[chord].first] = true;
                        has_chord[hor_chords_[chord].second] = true;
                    } else {
                        has_chord[ver_chords_[chord - hor_chords_.size()].first] = true;
                        has_chord[ver_chords_[chord - hor_chords_.size()].second] = true;
                    }
                
#ifdef DEBUG
                    if (chord < (int) hor_chords_.size()) {
                        std::cout << "Adding hor chord " << chord << " (" << hor_chords_[chord].first << ", " << hor_chords_[chord].second << ") == ";
                        std::cout << "(" << vertices_[hor_chords_[chord].first].first << ", " << vertices_[hor_chords_[chord].first].second << ") - (" << vertices_[hor_chords_[chord].second].first << ", " << vertices_[hor_chords_[chord].second].second << ") to final chord set " << std::endl; 
                    } else {
                        std::cout << "Adding ver chord " << chord << " (" << ver_chords_[chord - hor_chords_.size()].first << ", " << ver_chords_[chord - hor_chords_.size()].second << ") == ";
                        std::cout << "(" << vertices_[ver_chords_[chord - hor_chords_.size()].first].first << ", " << vertices_[ver_chords_[chord - hor_chords_.size()].first].second << ") - (" << vertices_[ver_chords_[chord - hor_chords_.size()].second].first << ", " << vertices_[ver_chords_[chord - hor_chords_.size()].second].second << ") to final chord set " << std::endl; 
                    }
#endif
                    break;
                }
            }
        }

        if (!has_chord[concave_vertices_id_[i]]) {
            std::cout << "\n!Found no viable chord for vertex " << concave_vertices_id_[i] << " (" << vertices_[concave_vertices_id_[i]].first << ", " << vertices_[concave_vertices_id_[i]].second << ")" << std::endl;
            success_ = false;
        }
    }

    for (int i = 0; i < (int) independent_chords_.size(); i++) {
        if (independent_chords_[i] < H) {
            final_chord_set.push_back(independent_chords_[i]);
        } else {
            final_chord_set.push_back(independent_chords_[i] + (hor_chords_.size() - H));
        }
    }    
    
    final_chord_set_ = final_chord_set;
}


Orientation RectilinearPartitioner::DetermineOrientation(int v1, int v2) {
    int hdiff, vdiff;
    hdiff = std::abs(vertices_[v2].first - vertices_[v1].first);
    vdiff = std::abs(vertices_[v2].second - vertices_[v1].second);

    if (hdiff == 0) {
        return HORIZONTAL;
    } else if (vdiff == 0) {
        return VERTICAL;
    } else {
        std::cout << "Something wrong in direction" << std::endl;
        success_ = false;
    }
    return HORIZONTAL;
}


std::vector<rect> RectilinearPartitioner::GetRectangleCoordinates() {
    // Given the chords of the partitioning, we need to find the rectangles

    // Chords can be part of at most two rectangles
    // Therefore pick a chord that has not been marked twice
    // starting point is first coordinate of chord, end point is second coordinate
    // While we have not found the end point
    //      Traverse the vertices where for every vertex 
    //      If it is convex, store it as corner
    //      Second we check if it has a chord, 
    //          If it does than we have to traverse it and mark it as visited

    // To achieve this the vertices are already in a useful order, we need to easily be able to know if they are convex, concave or rectilinear
    // Additionally, we need to know whether a chord is connected to it given that chord, we need to know the coordinates of its end points.

    std::vector<rect> rectangles;
    int nr_of_vertices = vertices_.size();
    int H = hor_chords_.size();
    int V = ver_chords_.size();
    int nr_of_chords = final_chord_set_.size();
    int chords_visited[H + V];  // -1 not used chord, 0 = not visited, 1 = visited anti-clockwise, 2 visited clockwise, 3 = visited clockwise and anti-clockwise
    memset(chords_visited, -1, sizeof(chords_visited));

    // set selected independent chords to 0
    for (int i = 0; i < (int) final_chord_set_.size(); i++) {
        chords_visited[final_chord_set_[i]] = 0;
    }

    int problem_count = 0;
    int chord_search_start = 0;
    while (true) {
        int direction = 1; // start anti-clockwise
        int start_vertex, end_vertex, current_vertex;
        start_vertex = end_vertex = -1;
        // find chord that has not been visited twice
        for (int i = chord_search_start; i < H + V; i++) {
            if (chords_visited[i] != -1 && chords_visited[i] < 3) {
                chord_search_start = i; // ensure next search starts from this vertex as all previous have been covered already
                // Get chord
                std::pair<int, int> chord;
                if (i < H) {
                    chord = hor_chords_[i];
                } else {
                    chord = ver_chords_[i - H];
                }

                if (chords_visited[i] == 0 || chords_visited[i] == 2) {
                    // default chord.first < chord.second
                    start_vertex = chord.first;
                    end_vertex = chord.second;
                } else {
                    // swap to ensure chord has been traversed in both ways
                    direction = 2; 
                    start_vertex = chord.second;
                    end_vertex = chord.first;
                }

                chords_visited[i] += direction;  // mark visited
                break;
            }
        }

        // if no chord found, we should be done
        if (start_vertex == -1) {
            break;
        }

        // vertex traversing always done by incrementing it
        current_vertex = (start_vertex+1)%nr_of_vertices;

        std::vector<std::pair<float, float>> corners;
        int previous_dir = DetermineOrientation(start_vertex, end_vertex);
        int previous_chord = -1;
        int previous_vertex = start_vertex;
        while (true) {
            if (problem_count > 10e5) {
                success_ = false;
                break;
            }
            // Every convex vertex is a corner
            int dir = DetermineOrientation(previous_vertex, current_vertex);
#ifdef DEBUG
            std::cout << previous_vertex << " - " << current_vertex << " prev_dir " << previous_dir << " cur_dir " << dir << std::endl;
#endif
            // Check if direction changed
            if (dir != previous_dir) {
                corners.push_back(std::make_pair(vertices_[current_vertex].first, vertices_[current_vertex].second));
            }
            previous_dir = dir;

            if (current_vertex == end_vertex) {
                break;
            }

            // if current_vertex part of other chord, follow this chord
            std::vector<int> vertex_chords = vertex_to_chord_[current_vertex];
            bool found_chord = false;
            for (int i = 0; i < (int) vertex_chords.size(); i++) {
                // If chord has not been visited in this direction yet
                if (chords_visited[vertex_chords[i]] != -1 && chords_visited[vertex_chords[i]] != 3 && vertex_chords[i] != previous_chord) {
                    // Extract chord
                    std::pair<int, int> chord;
                    if (vertex_chords[i] >= H) {
                        chord = ver_chords_[vertex_chords[i] - H];
                    } else {
                        chord = hor_chords_[vertex_chords[i]];
                    }

                    if (current_vertex == chord.first && chords_visited[vertex_chords[i]] != 2) {
                        chords_visited[vertex_chords[i]] += 2;  // mark traversed with second vertex first
                    } else if (current_vertex == chord.second && chords_visited[vertex_chords[i]] != 1) {
                        chords_visited[vertex_chords[i]] += 1;  // mark traversed in default direction
                    } else {
                        //std::cout << "Chord " << vertex_chords[i] << " rejected cause already traversed in specified direction" << std::endl;
                        problem_count++;
                    }

                    // set vertex for next iteration to other end of the chord
                    previous_vertex = current_vertex;
#ifdef DEBUG
                    std::cout << "Via chord" << std::endl;
#endif
                    if (current_vertex == chord.first) {
                        current_vertex = chord.second;
                    } else {
                        current_vertex = chord.first;
                    }

                    // indicate chord is found and 
                    found_chord = true;
                    previous_chord = vertex_chords[i];
                    break;
                }
                else {
#ifdef DEBUG
                    /*std::cout << "From current vertex " << current_vertex << " with chords ";
                    for (auto chord : vertex_chords) {
                        std::cout << chord << ", ";
                    }
                    std::cout << " no viable found" << std::endl;*/
#endif
                }
            }

            // else simply increment
            if (!found_chord) {
                previous_vertex = current_vertex;
                current_vertex = (current_vertex + 1) % nr_of_vertices;
            }
        }
#ifdef DEBUG
        std::cout << "#corners after loop " << corners.size() << std::endl;
#endif
        if (corners.size() == 4) {
            rectangles.push_back(CornersToRectangle(corners));
        } else if (corners.size() == 3) {
            corners.push_back(std::make_pair(vertices_[start_vertex].first, vertices_[start_vertex].second));
            rectangles.push_back(CornersToRectangle(corners));
        } else {
            std::cout << "!Error: not 3 or 4 corners rect (" << corners.size() << ")" << std::endl;
            success_ = false;
        }
    }


    /*std::cout << "vertices " << vertices_.size() << std::endl;
    std::cout << "convex vertices " << convex_vertices_id_.size() << std::endl;
    std::cout << "concave vertices " << concave_vertices_id_.size() << std::endl;
    std::cout << "straight vertices " << straight_vertices_id_.size() << std::endl;
    std::cout << "hor chords " << hor_chords_.size() << std::endl;
    std::cout << "ver chords " << ver_chords_.size() << std::endl;
    std::cout << "matching hor chords " << matching_hor_chords_.size() << std::endl;
    std::cout << "matching ver chords " << matching_ver_chords_.size() << std::endl;
    std::cout << "vertex to chord " << vertex_to_chord_.size() << std::endl;
    std::cout << "vertex to matched chord " << vertex_to_matched_chord_.size() << std::endl;
    std::cout << "independent chords " << independent_chords_.size() << std::endl;
    std::cout << "final chord set " << final_chord_set_.size() << std::endl;*/

    return rectangles;
}

rect RectilinearPartitioner::CornersToRectangle(std::vector<std::pair<float, float>> corners) {
    rect rectangle;

    auto xmin = *std::min_element(corners.cbegin(), corners.cend(), [](const auto& lhs, const auto& rhs) {
                        return lhs.first < rhs.first;    
                    });
    rectangle.x1 = xmin.first;

    auto ymin = *std::min_element(corners.cbegin(), corners.cend(), [](const auto& lhs, const auto& rhs) {
                        return lhs.second < rhs.second;    
                    });
    rectangle.y1 = ymin.second;

    auto xmax = *std::min_element(corners.cbegin(), corners.cend(), [](const auto& lhs, const auto& rhs) {
                        return lhs.first > rhs.first;    
                    });
    rectangle.x2 = xmax.first;

    auto ymax = *std::min_element(corners.cbegin(), corners.cend(), [](const auto& lhs, const auto& rhs) {
                        return lhs.second > rhs.second;    
                    });
    rectangle.y2 = ymax.second;

    return rectangle;
}
