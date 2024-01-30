#include "Statistics.h"

#include <iostream>


StatisticsWriter::StatisticsWriter(std::string file_name) {
    // Try to open file
    //file_ = std::ifstream(file_name, std::ios::out);
    file_.open(file_name, std::ios::out);
    if (!file_.is_open()) {
        std::cout << "Error creating statistics file: " << file_name << std::endl;
        exit(1);
    } else {
        std::cout << "Created statistics file: " << file_name << std::endl;
    }

    // Write header to file based on struct
    file_ << "grid_size,"
            << "selection_size,"
            << "shape,"
            << "nr_of_vertices,"
            << "nr_of_rectangles,"
            << "partition_time (ms),"
            << "query_method,"
            << "query_time (ms),"
            << "query_answer,"
            << "groundtruth,"
            << "avg_error_N,"
            << "avg_error_L1,"
            << "avg_error_L2,"
            << "avg_error_Rel,"
            << "epsilon,"
            << "memory,"
            << "memory_limit,"
            << "resolution,"
            << "FPR,"
            << "FNR,"
            << "precision,"
            << "recall,"
            << "f1_score,"
            << std::endl;
}

StatisticsWriter::~StatisticsWriter() {
    if (file_.is_open()) {
        file_.flush();
        file_.close();
    }
}

// Write statistics to file with comma as delimiter
void StatisticsWriter::WriteStatistics(statistics stats) {
    // Todo: there probably exists a generic method for this
    file_ << stats.grid_size << ","
            << stats.selection_size << ","
            << stats.shape + ","
            << stats.nr_of_vertices << ","
            << stats.nr_of_rectangles << ","
            << stats.partition_time << ","
            << stats.query_method << ","
            << stats.query_time << ","
            << stats.query_answer << ","
            << stats.groundtruth << ","
            << stats.avg_error_N << ","
            << stats.avg_error_L1 << ","
            << stats.avg_error_L2 << ","
            << stats.avg_error_Rel << ","
            << stats.epsilon << ","
            << stats.memory << ","
            << stats.memory_limit << ","
            << stats.resolution << ","
            << stats.false_positive_rate << ","
            << stats.false_negative_rate << ","
            << stats.precision << ","
            << stats.recall << ","
            << stats.f1_score << ","
            << std::endl;
}



UpdateStatisticsWriter::UpdateStatisticsWriter(std::string file_name) {
    
    // Try to open file
    file_ = std::ofstream(file_name, std::ios::out);
    //file_.open(file_name, std::ios::out);
    if (!file_.is_open()) {
        std::cout << "Error creating statistics file: " << file_name << std::endl;
        exit(1);
    } else {
        std::cout << "Created statistics file: " << file_name << std::endl;
    }
    file_ << "Ingestion time (ms),Memory (MB),Ingestion 100k times (ms)\n";
}

UpdateStatisticsWriter::~UpdateStatisticsWriter() {
    if (file_.is_open()) {
        file_.flush();
        file_.close();
    }
}

// Write statistics to file with comma as delimiter
void UpdateStatisticsWriter::WriteInsertTime(time_t time, uint memory, std::vector<time_t> sub_times) {
    // Todo: there probably exists a generic method for this
    file_ << time << "," << memory << ",[";
    for (int i = 0; i < (int) sub_times.size(); i++) {
        file_ << sub_times[i];
        if (i < (int) sub_times.size() - 1) file_ << ",";
    } 
    file_ << "]" << std::endl;
}

// Write statistics to file with comma as delimiter
void UpdateStatisticsWriter::WriteInsertTime(time_t time, std::string memory, std::vector<time_t> sub_times) {
    // Todo: there probably exists a generic method for this
    file_ << time << "," << memory << ",[";
    for (int i = 0; i < (int) sub_times.size(); i++) {
        file_ << sub_times[i];
        if (i < (int) sub_times.size() - 1) file_ << ",";
    } 
    file_ << "]" << std::endl;
}