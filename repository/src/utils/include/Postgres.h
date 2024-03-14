#ifndef POSTGRES_H_
#define POSTGRES_H_

#include <Utils.h>

#include <pqxx/pqxx>  // cpp postgresql api

typedef struct input_data {
    long timestamp;
    long ip;
    int x;
    int y;
} input_data;

typedef struct QuerySet {
    std::vector<range> ranges = {};
    std::vector<std::pair<float, float>> vertices = {};
    long item = 0;  // ip start
    long item_end = 0;  // ip end
    long groundtruth = -1;
    int L1 = -1;
    int N = -1;
    int offset_x = 0;
    int offset_y = 0;
    int time_window = 0;  // time window consting of [time_window, max_window], i.e., 0 is the entire time frame

    QuerySet() {};
    QuerySet(std::vector<range> ranges, std::vector<std::pair<float, float>> vertices, long item, long item_end, long groundtruth, int L1, int N, int offset_x, int offset_y, int time_window) : ranges(ranges), vertices(vertices), item(item), item_end(item_end), groundtruth(groundtruth), L1(L1), N(N), offset_x(offset_x), offset_y(offset_y), time_window(time_window) {};
} QuerySet;


class Postgres {
    public:
        // Datapath for bulk insertion
        Postgres(std::string data_path = "", int index = 1);
        ~Postgres();

        // Get size of entire relation + index(es)
        std::string GetDBSize();
        int GetNrRecords();

        // Insert single record
        void Insert(std::tuple<long, long, int, int> item, float &time);
        void InsertSet(std::vector<std::tuple<long, long, int, int>> data, time_t &time);
        void InsertSet(input_data* data, int count, time_t &time, std::vector<time_t> &interval_times);

        // COUNT(*)
        long RangeCount(std::vector<range> ranges, float &time);
        long RangeCount(std::vector<std::pair<float, float>> vertices, float &time);

        // SUM(packet_count)
        //int RangeSum(std::vector<range> ranges, time_t &time);
        //int RangeSum(std::vector<std::pair<float, float>> vertices, time_t &time);

        // Count(*) WHERE int_ip = item
        long RangeCountIP(std::vector<range> ranges, long item, float &time);
        long RangeCountIP(std::vector<std::pair<float, float>> vertices, long item, float &time);

        std::vector<long> QuerySetRangeCountIP(std::vector<QuerySet> queries, int index, float &time, int sample_size = 1);

        std::vector<long> GetItemsInRange(std::vector<range> ranges, int min_query_answer = 0);
        std::vector<QuerySet> GetIPRangeQueriesFrequency(std::vector<range> zerod_ranges, std::vector<std::pair<float, float>> vertices, int N, int resolution, int max_x_offset, int max_y_offset, int sample_size, bool range_queries, int min_query_answer=0, bool timestamp=false);
        std::vector<QuerySet> GetIPRangeQueriesCountDistinct(std::vector<range> zerod_ranges, std::vector<std::pair<float, float>> vertices, int N, int max_x_offset, int max_y_offset, int sample_size, int limit=-1, int min_query_answer=0);
        std::vector<QuerySet> GetIPRangeQueriesMembership(std::vector<range> zerod_ranges, std::vector<std::pair<float, float>> vertices, int N, int max_x_offset, int max_y_offset, int sample_size, int limit=-1, int min_query_answer=0);
        std::vector<QuerySet> GetIPRangeQueriesL2(std::vector<range> zerod_ranges, std::vector<std::pair<float, float>> vertices, int N, int max_x_offset, int max_y_offset, int sample_size, int min_query_answer=0);

    private:
        int index_ = 0;
        std::string ComposeRange(std::vector<range> ranges);
        std::string ComposeRange(std::vector<std::pair<float, float>> vertices);
        long QuerySingle(std::string query, float &time);
};

#endif  // POSTGRES_H_