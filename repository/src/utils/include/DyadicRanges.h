#ifndef DYADIC_RANGES_H_
#define DYADIC_RANGES_H_


typedef struct dyadic1D {
    long start = 0;
    long end = 0;
    float coverage = 1.0f;

    dyadic1D(long start = 0, long end = 0, float coverage = 1.0f) : start(start), end(end), coverage(coverage) {}
} dyadic1D;

typedef struct dyadic2D {
    int x1 = 0;
    int y1 = 0;
    int x2 = 0;
    int y2 = 0;
    float coverage = 1.0f;

    dyadic2D(int x1 = 0, int y1 = 0, int x2 = 0, int y2 = 0, float coverage = 1.0f) : x1(x1), y1(y1), x2(x2), y2(y2), coverage(coverage) {}
} dyadic2D;

typedef struct dyadic3D {
    int x1 = 0;
    int y1 = 0;
    int x2 = 0;
    int y2 = 0;
    long ip_start = 0;
    long ip_end = 0;
    float coverage = 1.0f;

    dyadic3D(int x1 = 0, int y1 = 0, int x2 = 0, int y2 = 0, long ip_start = 0, long ip_end = 0, float coverage = 1.0f) : x1(x1), y1(y1), x2(x2), y2(y2), ip_start(ip_start), ip_end(ip_end), coverage(coverage) {}
} dyadic3D;


#endif  // DYADIC_RANGES_H_