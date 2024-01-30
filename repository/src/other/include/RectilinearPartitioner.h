#ifndef RECTILINEAR_PARTITIONER_H_
#define RECTILINEAR_PARTITIONER_H_

#include "Utils.h"

#include <vector>

/* OLD RECTILINEAR POLYGON PARTITIONING CODE, CORRECT AND MORE EFFICIENT CODE IN PARTITIONER.CPP */

enum Direction {
    UP,
    DOWN,
    LEFT,
    RIGHT
};

enum Orientation {
    HORIZONTAL,
    VERTICAL
};

class RectilinearPartitioner {
    public:
        // Constructor takes vertices as input sorted in anti-clockwise order
        RectilinearPartitioner(std::vector<std::pair<float, float>> vertices);
        ~RectilinearPartitioner();

        // Compute rectangular partitioning based on input method
        std::vector<rect> ComputePartitioning(int method = 0);

        inline bool PartitioningSuccess() { return success_; }


    private:
        std::vector<std::pair<float, float>> vertices_;
        std::vector<int> convex_vertices_id_, concave_vertices_id_, straight_vertices_id_;
        std::vector<std::tuple<float, float, int>> concave_vertices_;
        std::vector<std::pair<int, int>> hor_chords_, ver_chords_, matching_hor_chords_, matching_ver_chords_;
        std::vector<std::vector<int>> vertex_to_chord_, vertex_to_matched_chord_;
        std::vector<int> independent_chords_, unusable_chords_, final_chord_set_;
        bool success_ = true;

        // Compute the rectangular partitioning based on the paper by Imai and Asano, O(V^(3/2)log(n))
        std::vector<rect> ComputePartitioningMNC();

        // Compute the rectangular partitioning based on taking the minimum number of horizontal or vertical chords.
        // Note, only produces correct output for orthogonally convex polygons
        std::vector<rect> ComputePartitioningMinChord();

        // Given that graph has only four convex vertices, return this as rectangle
        std::vector<rect> ReturnRect();

        // Compute whether vertex is convex, concave or straight
        bool ComputeVertexType();  

        // Compute horizontal and vertical chords in graph
        void ComputeChords();

        // Given chords, compute maximal non-intersecting chord set
        void MaximumIndependentSet();

        // Bipartite matching
        bool BipartiteMatch(bool **graph, int h, bool seen[], int matchH[], int matchV[]);

        // Given matching, compute chosen chords
        void FinalizeChordSet();

        // Given final chord set, compute the rectangles they create
        std::vector<rect> GetRectangleCoordinates();

        // Convert corners to rectangle range
        rect CornersToRectangle(std::vector<std::pair<float, float>> corners);

        // Given v1 and v2, determine whether direction is horiontal or vertical
        Orientation DetermineOrientation(int v1, int v2);

        // Given two chords, determine if they intersect
        bool CheckIntersect(std::pair<int, int> chord1, std::pair<int, int> chord2);


};

#endif  // RECTILINEAR_PARTITIONER_H_