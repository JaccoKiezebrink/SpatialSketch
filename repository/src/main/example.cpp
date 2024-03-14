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

    for (int i = 0; i < 1001; i++) {
        sp.Update(0, 0, ip, i);
        if (i % 100 == 0) {
            std::cout << i << " Memory usage: " << sp.GetSize() << std::endl;
        }
    }

    /*ECM ecm = ECM(0.1, 0.05, nullptr);
    ecm.Insert(1, 1);
    ecm.Insert(1, 2);
    ecm.Insert(1, 3);
    ecm.Insert(1, 4);
    ecm.Insert(1, 5);
    ecm.Insert(1, 6);

    std::cout << "ECM query answer to item 1 in window [6, 7] is " << ecm.QueryItem(1, 7) << std::endl;
    std::cout << "ECM query answer to item 1 in window [6, 6] is " << ecm.QueryItem(1, 6) << std::endl;
    std::cout << "ECM query answer to item 1 in window [5, 6] is " << ecm.QueryItem(1, 5) << std::endl;
    std::cout << "ECM query answer to item 1 in window [4, 6] is " << ecm.QueryItem(1, 4) << std::endl;
    std::cout << "ECM query answer to item 1 in window [3, 6] is " << ecm.QueryItem(1, 3) << std::endl;
    std::cout << "ECM query answer to item 1 in window [2, 6] is " << ecm.QueryItem(1, 2) << std::endl;
    std::cout << "ECM query answer to item 1 in window [1, 6] is " << ecm.QueryItem(1, 1) << std::endl;
    std::cout << "ECM query answer to item 1 in window [0, 6] is " << ecm.QueryItem(1, 0) << std::endl;*/


    return 0;
}