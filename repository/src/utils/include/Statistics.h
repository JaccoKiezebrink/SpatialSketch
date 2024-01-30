#ifndef STATISTICS_H_
#define STATISTICS_H_

#include <string>
#include <fstream>
#include <chrono>
#include <vector>

// Struct of statistics to be logged
typedef struct statistics {
    // Basics
    int grid_size = -1;
    int selection_size = -1;
    std::string shape;
    int nr_of_vertices = -1;

    // Rectilinear partitioning
    int nr_of_rectangles = -1;
    float partition_time = -1;

    // Query
    std::string query_method;
    int nr_of_subqueries = -1;
    float query_time = -1;
    long query_answer = -1;
    long groundtruth = -1;
    float avg_error_N = -1;
    float avg_error_L1 = -1;
    float avg_error_L2 = -1;
    float avg_error_Rel = -1;

    // Membership queries
    float false_positive_rate = -1;
    float false_negative_rate = -1;
    float precision = -1;
    float recall = -1;
    float f1_score = -1;

    // Stats
    float epsilon = -1;
    std::string memory;
    std::string memory_limit;
    int resolution = -1;


} statistics;


class StatisticsWriter {
    public:
        StatisticsWriter(std::string file_name);
        ~StatisticsWriter();
        void WriteStatistics(statistics stats);
    private:
        std::ofstream file_;
};

class UpdateStatisticsWriter {
    public:
        UpdateStatisticsWriter(std::string file_name);
        ~UpdateStatisticsWriter();
        void WriteInsertTime(time_t time, uint memory = 0, std::vector<time_t> sub_times = std::vector<time_t>());
        void WriteInsertTime(time_t time, std::string memory, std::vector<time_t> sub_times = std::vector<time_t>());
    private:
        std::ofstream file_;
};

#endif  // STATISTICS_H_