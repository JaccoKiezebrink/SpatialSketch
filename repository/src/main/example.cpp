#include "SpatialSketch.h"

#include <iostream>

int main() {
    // Initialize SpatialSketch with desired sketch, grid size, and optionally a memory limit, epsilon, delta, and domain size
    SpatialSketch sp = SpatialSketch("CM", 1024);

    // Insert some ip address at given x, y
    long ip = 101;
    sp.Update(1, 1, ip, 1);
    sp.Update(2, 2, ip, 3);

    // Determine range to query
    range r = range(0, 0, 3, 3);
    std::vector<range> ranges = {r};  // SpatialSketch will handle multiple orthogonal ranges if given
    long query_answer = sp.QueryRanges(ranges, ip);

    std::cout << "\nQuery answer to range " << r.x1 << ", " << r.y1 << ", " << r.x2 << ", " << r.y2 << " is " << query_answer << std::endl;

    return 0;
}