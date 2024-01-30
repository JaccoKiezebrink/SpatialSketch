#ifndef PARTITIONER_H_
#define PARTITIONER_H_

#include "Utils.h"

#include <vector>
#include <set>


enum Direction {
    DIR_UNKNOWN = 0,
    UP,
    DOWN,
    LEFT,
    RIGHT
};

enum VertexType {
    TYPE_UNKNOWN = 0,
    CONVEX,
    CONCAVE,
    STRAIGHT
};

typedef struct vertex {
    int id;
    float x;
    float y;
    VertexType type;
    Direction dir;

    vertex(int id, float x, float y) : id(id), x(x), y(y) {}
    vertex(int id, float x, float y, VertexType type) : id(id), x(x), y(y), type(type) {}
    vertex(int id, float x, float y, VertexType type, Direction dir) : id(id), x(x), y(y), type(type), dir(dir) {}
} vertex;


class Partitioner {
    public:
        // Constructor takes vertices as input sorted in anti-clockwise order
        Partitioner(std::vector<std::pair<float, float>> vertices);
        //~RectilinearPartitioner();

        // Compute rectangular partitioning
        std::vector<rect> ComputePartitioning();

        // Indicate whether any error occured during partitioning
        inline bool PartitioningSuccess() {
            return partition_success_;
        }

        void PrintRectangles(std::vector<rect> rectangles);
        
    private:
        int N_;
        int concave_count_;
        bool partition_success_ = true;
        std::vector<std::pair<float, float>> orig_vertices_;
        std::vector<vertex> vertices_;
        std::vector<std::pair<vertex, vertex>> min_chords_;
        std::vector<std::pair<int, int>> min_chords_id_;

        std::set<std::pair<std::pair<float, float>, std::pair<float, float>>> vertex_map_;

        Direction DetermineDirection(std::pair<float, float> v, std::pair<float, float> w);
        Direction DetermineDirection(vertex v, vertex w);

        // Compute the vertex types and seperate vectors per type
        bool ComputeVertexType();

        // Determine whether two vertices lay adjacent to each other
        bool AreNeighbors(vertex v, vertex w);

        // Compute chords and return smallest of horizontal or vertical set
        bool ComputeMinChords();

        std::vector<vertex> ReduceVertices(std::vector<vertex> vertices);
        std::vector<rect> CornerIdentification(std::vector<vertex> vertices);
};

#endif  // PARTITIONER_H_