#include "SpatialSketch.h"
#include "ECM.h"

#include <iostream>

int main() {
    // Initialize SpatialSketch with desired sketch, grid size, and optionally a memory limit, epsilon, delta, and domain size
    SpatialSketch sp = SpatialSketch("ECM", 1024, 3350000);
    std::cout << "Memory usage: " << sp.GetSize() << std::endl;

    // Insert some ip address at given x, y
    long ip = 101;
    sp.Update(0, 0, ip, 0);
    sp.Update(0, 0, ip, 1);
    sp.Update(0, 0, ip, 2);
    sp.Update(0, 0, ip, 3);
    sp.Update(0, 0, ip, 4);

    // Determine range to query
    range r = range(0, 0, 3, 3);
    std::vector<range> ranges = {r};  // SpatialSketch will handle multiple orthogonal ranges if given
    long query_answer = sp.QueryRanges(ranges, ip, ip, 0);

    std::cout << "\nQuery answer to range " << r.x1 << ", " << r.y1 << ", " << r.x2 << ", " << r.y2 << " is " << query_answer << std::endl;

    return 0;
}