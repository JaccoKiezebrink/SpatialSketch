#include "Postgres.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <random>

//#define PRINT_QUERIES
#define KUBERNETES  // Define used when running on kubernetes cluster which changes bulk insert from file to stdin

Postgres::Postgres(std::string data_path, int index) {
    index_ = index;

    // Try to open psql connection
    pqxx::connection C("dbname = postgres user = postgres password = postgres hostaddr = 127.0.0.1 port = 5432");
    if (C.is_open()) {
        // Get time before creating table and inserting data
        timespec timesp;
        time_t before_time;
        clock_gettime(CLOCK_REALTIME, &timesp);
        before_time = timesp.tv_sec * 1000 + timesp.tv_nsec / 1000000;


        pqxx::work w(C);
        // Cancel all active queries
        w.exec("SELECT pg_cancel_backend(pid) FROM pg_stat_activity WHERE state = 'active' and pid <> pg_backend_pid();");
        // Drop existing tables and views
        w.exec("DROP MATERIALIZED VIEW IF EXISTS grouped_ip;");
        w.exec("DROP MATERIALIZED VIEW IF EXISTS grouped_ip_cd;");
        w.exec("DROP MATERIALIZED VIEW IF EXISTS grouped_ip_mem;");
        w.exec("DROP MATERIALIZED VIEW IF EXISTS distinct_ips;");


        w.exec("DROP TABLE IF EXISTS iplocations;");

        // Rebuild the table, insert data, analyze and build any index if specified
        w.exec("CREATE UNLOGGED TABLE iplocations (timestamp bigint, ip bigint, longitude int, latitude int);");
        // Bulk insert data if data_path is not empty
        if (data_path != "") {
#ifdef KUBERNETES
            std::cout << "Postgres STDIN copy start" << std::endl;
            std::ifstream* data_file = new std::ifstream(data_path, std::ios::in);
            if (!data_file->is_open()) {
                std::cout << "Error opening data file: " << data_path << std::endl;
                return;
            }

            //w.exec("COPY iplocations FROM STDIN WITH DELIMITER ',';");
            pqxx::stream_to S(w, "iplocations");
            std::string line, cell;
            std::getline(*data_file,line); // header
            long timestamp = 0;
            while (std::getline(*data_file, line)) {
                timestamp++;
                std::stringstream linestream(line);
                std::tuple<long, long, int, int> tup;
                for (int i = 0; i < 4; i++) {
                    std::getline(linestream, cell, ',');
                    // TODO: This depends on the order of the columns in the dataset, therefore not ideal
                    // Current GeoCaide orders it as follows [timestamp, int_ip, long, lat, points]
                    if (i == 0) std::get<0>(tup) = timestamp; //(long) std::stol(cell); // timestamp
                    if (i == 1) std::get<1>(tup) = (long) std::stol(cell);  // int_ip (can be 2^32, thus uint/long)
                                                                            // We cast it to int as current fucntion support this, that does imply negative items can occur, but this should not change any result.
                    if (i == 2) std::get<2>(tup) = std::stoi(cell);  // long
                    if (i == 3) std::get<3>(tup) = std::stoi(cell);  // lat
                }
                S << tup;
            }
            S.complete();
            data_file->close();

#else
            // Insert data
            std::string insert_command = "COPY iplocations(timestamp, ip, longitude, latitude) FROM '" + data_path + "' DELIMITER ',' CSV HEADER;";
            w.exec(insert_command);
#endif
        }
        
        // Analyze inserted data
        w.exec("ANALYZE iplocations;");

        // Build index based on parsed parameter
        if (index == 1) {
            w.exec("CREATE INDEX sorted_ind ON iplocations USING btree (longitude, latitude);");
        } else if (index == 2) {
            w.exec("CREATE INDEX sorted_ind2 ON iplocations USING btree (ip, longitude, latitude);");
        } else if (index == 3 || index == 4) {
            w.exec("CREATE INDEX gist_ind ON iplocations USING gist (point(longitude, latitude));");
        } else if (index == 5 || index == 6) {
            w.exec("CREATE INDEX spgist_ind ON iplocations USING spgist (point(longitude, latitude));");
        }
        w.commit();

        // Log end time
        clock_gettime(CLOCK_REALTIME, &timesp);
        std::cout << "Opened Postgres database and inserted data in " << ((timesp.tv_sec * 1000 + timesp.tv_nsec / 1000000) - before_time) << "ms" << std::endl;
    
        pqxx::work qmem(C);
        pqxx::result result = qmem.exec("SELECT pg_size_pretty(pg_total_relation_size('iplocations'));");
        std::cout << "Postgres table + index size: " << result[0][0].c_str() << std::endl << std::endl;
    } else {
        std::cout << "Could not open Postgres connection" << std::endl;
    }
}

Postgres::~Postgres() {
    // Drop table
    try {
        pqxx::connection C("dbname = postgres user = postgres password = postgres hostaddr = 127.0.0.1 port = 5432");
        if (C.is_open()) {
            try {
                pqxx::work dropmview(C);
                dropmview.exec("DROP MATERIALIZED VIEW IF EXISTS grouped_ip;");
                dropmview.exec("DROP MATERIALIZED VIEW IF EXISTS grouped_ip_cd;");
                dropmview.exec("DROP MATERIALIZED VIEW IF EXISTS grouped_ip_mem;");
                dropmview.exec("DROP MATERIALIZED VIEW IF EXISTS distinct_ips;");

                dropmview.commit();

                pqxx::work w(C);
                w.exec("DROP TABLE IF EXISTS iplocations;");
                w.commit();
            } catch (const std::exception &e) {
                std::cerr << "Exception at Postgres deconstructor: " << e.what() << std::endl;
            }
            C.disconnect();
        }
    } catch (const std::exception &e) {
        std::cerr << "Exception at Postgres desconstructor:" << e.what() << std::endl;
    }
}

std::string Postgres::GetDBSize() {
    try {
        // Todo: Possibly it is not required to open connection each and every time
        pqxx::connection C("dbname = postgres user = postgres password = postgres hostaddr = 127.0.0.1 port = 5432");
        if (C.is_open()) {
            pqxx::work qmem(C);
            pqxx::result result = qmem.exec("SELECT pg_size_pretty(pg_total_relation_size('iplocations'));");
            //std::cout << "Postgres table + index size: " << result[0][0].c_str() << std::endl << std::endl; 
            return result[0][0].c_str();
        }
    } catch (const std::exception &e) {
        std::cerr << "Exception at Insert: " << e.what() << std::endl;
    }
    return "unknown";
}

int Postgres::GetNrRecords() {
    try {
        // Todo: Possibly it is not required to open connection each and every time
        pqxx::connection C("dbname = postgres user = postgres password = postgres hostaddr = 127.0.0.1 port = 5432");
        if (C.is_open()) {
            pqxx::work qmem(C);
            pqxx::result result = qmem.exec("SELECT COUNT(*) from iplocations;");
            //std::cout << "Postgres table + index size: " << result[0][0].c_str() << std::endl << std::endl; 
            return result.at(0)[0].as<int>();
        }
    } catch (const std::exception &e) {
        std::cerr << "Exception at Insert: " << e.what() << std::endl;
    }
    return -1;
}


void Postgres::Insert(std::tuple<long, long, int, int> item, float &time) {
    try {
        // Todo: Possibly it is not required to open connection each and every time
        pqxx::connection C("dbname = postgres user = postgres password = postgres hostaddr = 127.0.0.1 port = 5432");
        if (C.is_open()) {
            pqxx::work w(C);

            std::string query = "INSERT INTO iplocations (timestamp, ip, longitude, latitude) VALUES ('" + std::to_string(std::get<0>(item)) + "', '" + std::to_string(std::get<1>(item)) + "', '" + std::to_string(std::get<2>(item)) + "', '" + std::to_string(std::get<3>(item)) + "');";

            // Note that time is measured within connect block to avoid measuring connection time
            timespec timesp;
            time_t before_time;
            clock_gettime(CLOCK_REALTIME, &timesp);
            before_time = timesp.tv_sec * 1000000000 + timesp.tv_nsec;
            
            // Execute query
            w.exec(query);

            // end time
            clock_gettime(CLOCK_REALTIME, &timesp);
            time = ((timesp.tv_sec * 1000000000 + timesp.tv_nsec) - before_time) / 1000;     
        }
    } catch (const std::exception &e) {
        std::cerr << "Exception at Insert: " << e.what() << std::endl;
    }
}

void Postgres::InsertSet(std::vector<std::tuple<long, long, int, int>> data, time_t &time) {
    timespec before_time, after_time;
    int insert_count = 0;
    std::string query;
    try {
        // Todo: Possibly it is not required to open connection each and every time
        pqxx::connection C("dbname = postgres user = postgres password = postgres hostaddr = 127.0.0.1 port = 5432");
        if (C.is_open()) {
            pqxx::work *w;
            w = new pqxx::work(C);

            // Note that time is measured within connect block to avoid measuring connection time
            
            clock_gettime(CLOCK_REALTIME, &before_time);

            for (auto tup : data) {
                query = "INSERT INTO iplocations (timestamp, ip, longitude, latitude) VALUES ('" + std::to_string(std::get<0>(tup)) + "', '" + std::to_string(std::get<1>(tup)) + "', '" + std::to_string(std::get<2>(tup)) + "', '" + std::to_string(std::get<3>(tup)) + "');";
                
                // Execute query
                w->exec(query);

                // Commit every 10k records
                insert_count++;
                if (insert_count % 10000 == 0) {
                    w->commit();
                    delete w;
                    w = new pqxx::work(C);
                }
            }
        }
        // end time
        C.disconnect();
        clock_gettime(CLOCK_REALTIME, &after_time);
        time = (after_time.tv_sec - before_time.tv_sec) * 1000 - (after_time.tv_nsec - before_time.tv_nsec) / 1000000;
    } catch (const std::exception &e) {
        std::cerr << "Exception at Insert: " << e.what() << std::endl;
    }
}


void Postgres::InsertSet(input_data* data, int count, time_t &time, std::vector<time_t> &interval_times) {
    timespec before_time, after_time;
    int insert_count = 0;
    std::string query;
    try {
        // Todo: Possibly it is not required to open connection each and every time
        pqxx::connection C("dbname = postgres user = postgres password = postgres hostaddr = 127.0.0.1 port = 5432");
        if (C.is_open()) {
            pqxx::work *w;
            w = new pqxx::work(C);

            // Note that time is measured within connect block to avoid measuring connection time
            time_t before_subtime;
            timespec subtimesp;
            clock_gettime(CLOCK_REALTIME, &before_time);
            clock_gettime(CLOCK_REALTIME, &subtimesp);
            before_subtime = subtimesp.tv_sec * 1000000 + subtimesp.tv_nsec / 1000;

            for (int i = 0; i < count; i++) {
                query = "INSERT INTO iplocations (timestamp, ip, longitude, latitude) VALUES ('" + std::to_string(data[i].timestamp) + "', '" + std::to_string(data[i].ip) + "', '" + std::to_string(data[i].x) + "', '" + std::to_string(data[i].y) + "');";
                
                // Execute query
                w->exec(query);

                // Commit every 10k records
                insert_count++;
                if (insert_count % 10000 == 0) {
                    w->commit();
                    delete w;
                    w = new pqxx::work(C);
                }
                 if ((insert_count) % 100000 == 0) {
                    clock_gettime(CLOCK_REALTIME, &subtimesp);
                    interval_times.push_back(((subtimesp.tv_sec * 1000000 + subtimesp.tv_nsec / 1000) - before_subtime) / 1000);
                    before_subtime = subtimesp.tv_sec * 1000000 + subtimesp.tv_nsec / 1000;
                }
            }
        }
        // end time
        C.disconnect();
        clock_gettime(CLOCK_REALTIME, &after_time);
        time = (after_time.tv_sec - before_time.tv_sec) * 1000 - (after_time.tv_nsec - before_time.tv_nsec) / 1000000;
    } catch (const std::exception &e) {
        std::cerr << "Exception at Insert: " << e.what() << std::endl;
    }
}

// Given set of ranges (rectangles) compose SQL range query to be appended to .. WHERE
std::string Postgres::ComposeRange(std::vector<range> ranges) {
    std::string range = "";

    if (index_ >= 3) {  // box
        for (size_t i = 0; i < ranges.size(); i++) {
            range += "point(longitude, latitude) <@ box '(" + std::to_string(ranges[i].x1) + "," + std::to_string(ranges[i].y1) + "),(" + std::to_string(ranges[i].x2) + "," + std::to_string(ranges[i].y2) + ")'";            
            if (i < ranges.size() - 1) {
                range += " OR ";
            }
        }
    }
    // Default, No index or btree
    else {
        for (size_t i = 0; i < ranges.size(); i++) {
            range += "(longitude >= " + std::to_string(ranges[i].x1) + " AND longitude <= " + std::to_string(ranges[i].x2) + " AND latitude >= " + std::to_string(ranges[i].y1) + " AND latitude <= " + std::to_string(ranges[i].y2) + ")";
            if (i < ranges.size() - 1) {
                range += " OR ";
            }
        }
    } 

    return range;
}

// Given polygon vertices compose SQL query to be appended to .. WHERE
std::string Postgres::ComposeRange(std::vector<std::pair<float, float>> vertices) {
    std::string polygon = "point(longitude, latitude) <@ polygon '(";
    for (size_t i = 0; i < vertices.size(); i++) {
        polygon += "(" + std::to_string(vertices[i].first) + "," + std::to_string(vertices[i].second) + "),";            
    }
    polygon += "(" + std::to_string(vertices[0].first) + "," + std::to_string(vertices[0].second) + "))'";
    return polygon;
}

// Given SQL query, execute it, return result and measure time
long Postgres::QuerySingle(std::string query, float &time) {
    try {
        // Todo: Possibly it is not required to open connection each and every time
        pqxx::connection C("dbname = postgres user = postgres password = postgres hostaddr = 127.0.0.1 port = 5432");
        if (C.is_open()) {
            pqxx::work w(C);
            pqxx::result result;

            // Explain analyze for timing
            /*result = w.exec("EXPLAIN ANALYZE " + query);
            
            // Loop over every row, if 'Time:' is in the string extract the runtime and accumulate
            float runtime = 0;
            for (auto const &row: result) {
                for (auto const &field: row) {
                    std::string sfield = field.c_str();
                    auto ptr = sfield.find("Time: ");
                    if (ptr != std::string::npos) {
                        sfield.erase(0, ptr + 6);  // Erase everything up to and including 'Time: '
                        runtime += std::stof(sfield.substr(0, sfield.find("ms") - 1));
                    }
                }
            }*/

            // Note that time is measured within connect block to avoid measuring connection time
            timespec timesp;
            time_t before_time;
            clock_gettime(CLOCK_REALTIME, &timesp);
            before_time = timesp.tv_sec * 1000000000 + timesp.tv_nsec;
            
            // Execute query
            result = w.exec(query);
            w.commit();

            // end time
            clock_gettime(CLOCK_REALTIME, &timesp);
            time = ((timesp.tv_sec * 1000000000 + timesp.tv_nsec) - before_time) / 1000000.f;     
            //std::cout << "Postgres time " << runtime << " ms vs. our time " << time << " ms" << std::endl; 

            // Get query result
            long res = 0;
            if (result.at(0)[0].size() != 0) {  // nan check
                res = result.at(0)[0].as<long>();
            }
            return res;
        }
    } catch (const std::exception &e) {
        std::cerr << "Exception at RangeCount: " << e.what() << std::endl;
    }
    return -1;  // return -1 if something went wrong
}

// Range count queries for ranges
long Postgres::RangeCount(std::vector<range> ranges, float &time) {
    std::string query = "SELECT COUNT(*) FROM iplocations WHERE " + ComposeRange(ranges) + ";";
#ifdef PRINT_QUERIES
    std::cout << "Query: " << query << std::endl;
#endif
    return QuerySingle(query, time);
}

// Range count queries for polygons
long Postgres::RangeCount(std::vector<std::pair<float, float>> vertices, float &time) {
    std::string query = "SELECT COUNT(*) FROM iplocations WHERE " + ComposeRange(vertices) + ";";
#ifdef PRINT_QUERIES
    std::cout << "Query: " << query << std::endl;
#endif
    return QuerySingle(query, time);
}

// Range sum queries for ranges
/*int Postgres::RangeSum(std::vector<range> ranges, time_t &time) {
    std::string query = "SELECT SUM(packet_count) FROM iplocations WHERE " + ComposeRange(ranges) + ";";
#ifdef PRINT_QUERIES
    std::cout << "Query: " << query << std::endl;
#endif
    return QuerySingle(query, time);
}*/

// Range count queries for polygons
/*int Postgres::RangeSum(std::vector<std::pair<float, float>> vertices, time_t &time) {
    std::string query = "SELECT SUM(packet_count) FROM iplocations WHERE " + ComposeRange(vertices) + ";";
#ifdef PRINT_QUERIES
    std::cout << "Query: " << query << std::endl;
#endif
    return QuerySingle(query, time);
}*/


long Postgres::RangeCountIP(std::vector<range> ranges, long item, float &time) {
    std::string query = "SELECT COUNT(*) FROM iplocations WHERE (" + ComposeRange(ranges) + ") AND ip = '" + std::to_string(item) + "';";
#ifdef PRINT_QUERIES
    std::cout << "Query: " << query << std::endl;
#endif
    return QuerySingle(query, time);
}

long Postgres::RangeCountIP(std::vector<std::pair<float, float>> vertices, long item, float &time) {
    std::string query = "SELECT COUNT(*) FROM iplocations WHERE (" + ComposeRange(vertices) + ") AND ip = '" + std::to_string(item) + "';";
#ifdef PRINT_QUERIES
    std::cout << "Query: " << query << std::endl;
#endif
    return QuerySingle(query, time);
}


std::vector<long> Postgres::QuerySetRangeCountIP(std::vector<QuerySet> queries, int index, float &time, int sample_size) {
    std::vector<long> answers = {};
    answers.reserve(queries.size());

    try {
        // Todo: Possibly it is not required to open connection each and every time
        pqxx::connection C("dbname = postgres user = postgres password = postgres hostaddr = 127.0.0.1 port = 5432");
        if (C.is_open()) {
            std::string query;

            // Note that time is measured within connect block to avoid measuring connection time
            timespec timesp;
            time_t before_time;
            clock_gettime(CLOCK_REALTIME, &timesp);
            before_time = timesp.tv_sec * 1000000000 + timesp.tv_nsec;

            for (auto q : queries) {
                for (int j = 0; j < sample_size; j++) {
                    if (index == 4 || index == 6) {
                        query = "SELECT COUNT(*) FROM iplocations WHERE (" + ComposeRange(q.vertices) + ") AND ip = '" + std::to_string(q.item) + "';";
                    } else {
                        query = "SELECT COUNT(*) FROM iplocations WHERE (" + ComposeRange(q.ranges) + ") AND ip = '" + std::to_string(q.item) + "';";
                    }

                    // Execute query
                    pqxx::work w(C);
                    pqxx::result result;
                    result = w.exec(query);
                    w.commit();

                    if (j == 0) {
                        answers.push_back(result.at(0)[0].as<long>());
                        //std::cout << result.at(0)[0].as<int>() << ", ";
                    }
                }
            }

            // obtain average time (milliseoncds) per query, by dividing total time by #samples and #queries
            clock_gettime(CLOCK_REALTIME, &timesp);
            time = ((timesp.tv_sec * 1000000000 + timesp.tv_nsec) - before_time) / (float) sample_size / (float) queries.size() / 1000000.f;    
            std::cout << std::endl; 
            C.disconnect();
        }
    } catch (const std::exception &e) {
        std::cerr << "Exception at RangeCount: " << e.what() << std::endl;
    }
    return answers;
}



std::vector<long> Postgres::GetItemsInRange(std::vector<range> ranges, int min_answer) {
    std::string query = "SELECT ip FROM iplocations WHERE (" + ComposeRange(ranges) + ") GROUP BY ip HAVING COUNT(ip) >= " + std::to_string(min_answer) + " ORDER BY COUNT(ip) DESC LIMIT 1;";
    try {
        // Todo: Possibly it is not required to open connection each and every time
        pqxx::connection C("dbname = postgres user = postgres password = postgres hostaddr = 127.0.0.1 port = 5432");
        if (C.is_open()) {
            pqxx::work w(C);
            pqxx::result result = w.exec(query);

            // Get query result
            std::vector<long> res;
            for (auto const &row: result) {
                for (auto const &field: row) {
                    res.push_back(std::stoul(field.c_str()));
                }
            }

            C.disconnect();
            return res;
        }
    } catch (const std::exception &e) {
        std::cerr << "Exception at GetIPInRange: " << e.what() << std::endl;
    }
    return {};  // return empty if something went wrong
}


long ScaleToResolutionCeil(long value, long current_resolution, long new_resolution) {
    return (long) ceil(new_resolution * (value / (float) current_resolution));
}

long ScaleToResolutionFloor(long value, long current_resolution, long new_resolution) {
    return (long) floor(new_resolution * (value / (float) current_resolution));
}


std::vector<QuerySet> Postgres::GetIPRangeQueriesFrequency(std::vector<range> zerod_ranges, std::vector<std::pair<float, float>> vertices, int N, int resolution, int max_x_offset, int max_y_offset, int sample_size, bool range_queries, int min_query_answer, bool timestamp) {
    std::vector<range> offset_ranges;
    std::vector<std::pair<float, float>> offset_vertices;
    std::vector<QuerySet> query_sets = {};
    query_sets.reserve(sample_size);
    int time_window = 0;
    pqxx::result neighbour_ips_result;
    try {
        // Todo: Possibly it is not required to open connection each and every time
        pqxx::connection C("dbname = postgres user = postgres password = postgres hostaddr = 127.0.0.1 port = 5432");
        if (C.is_open()) {
            // Get dataset size
            pqxx::work Nquery(C);
            int dataset_size = Nquery.exec("SELECT COUNT(*) FROM iplocations;")[0][0].as<int>();
            Nquery.commit();

            // Get top N items, their counts and location
            // By creating materialized view first
            pqxx::work mview(C);
            mview.exec("CREATE MATERIALIZED VIEW IF NOT EXISTS grouped_ip AS (SELECT ip, longitude, latitude, count(ip), max(timestamp) as timestamp FROM iplocations GROUP BY ip, longitude, latitude ORDER BY ip, longitude, latitude);");
            mview.exec("CREATE INDEX IF NOT EXISTS mview_btree ON grouped_ip USING btree (ip, longitude, latitude);");
            if (range_queries) {
                mview.exec("CREATE INDEX IF NOT EXISTS mview_btree2 ON grouped_ip USING btree (longitude, latitude);");            
            }
            if (timestamp) {
                mview.exec("CREATE INDEX IF NOT EXISTS btree_time ON iplocations USING btree (timestamp, longitude, latitude);");
            }
            mview.commit();
            // std::string q = "SELECT ip, count(ip) FROM grouped_ip GROUP BY ip HAVING count(ip) > 2;"; // opt 
            // pqxx::work w1(C);
            // pqxx::result res = w1.exec(q);
            // w1.commit();
            // for (int i = 0; i < dataset_size; i++) {
            //     std::cout << res[i][0].as<uint>() << std::endl;
            // }

            // Get all ip address that occur at atleast two locations
            std::string query = "SELECT ip, count(ip) FROM grouped_ip GROUP BY ip HAVING count(ip) > 2 ORDER BY count(ip) DESC LIMIT " + std::to_string(sample_size * 10) + ";"; // opt 
            pqxx::work w(C);
            pqxx::result result = w.exec(query);
            w.commit();

            int x_dim = (N - max_x_offset) - 1;
            int y_dim = (N - max_y_offset) - 1;
            int cmax_x_offset = max_x_offset; //ScaleToResolutionFloor(ScaleToResolutionFloor(N - x_dim, N, resolution), resolution, N);
            int cmax_y_offset = max_y_offset; //ScaleToResolutionFloor(ScaleToResolutionFloor(N - y_dim, N, resolution), resolution, N);

            int retries = 0;
            int samples = 0;
            int i = -1;
            while (samples < sample_size) {
                i++;
                pqxx::row entry;
                
                entry = result[(int) std::rand() % result.size()];
                
                long item = entry[0].as<long>();

                // Order differently to reduce bias towards some geolocation
                if (i % 4 == 0) {
                    query = "SELECT longitude, latitude FROM grouped_ip WHERE ip = " + std::to_string(item) + " ORDER BY longitude, latitude DESC;";
                } else if (i % 4 == 1) {
                    query = "SELECT longitude, latitude FROM grouped_ip WHERE ip = " + std::to_string(item) + " ORDER BY latitude, longitude DESC;";
                } else if (i % 4 == 2) {
                    query = "SELECT longitude, latitude FROM grouped_ip WHERE ip = " + std::to_string(item) + " ORDER BY longitude, latitude ASC;";
                } else {
                    query = "SELECT longitude, latitude FROM grouped_ip WHERE ip = " + std::to_string(item) + " ORDER BY latitude, longitude ASC;";
                }
                pqxx::work ips(C);
                pqxx::result ip_result = ips.exec(query);
                ips.commit();

                int x1, y1, x2, y2, x_dist, y_dist;
                int long_centroid, lat_centroid;
                for (int j = 1; j < (int) ip_result.size(); j++) {
                    x1 = ip_result[j-1][0].as<int>();
                    y1 = ip_result[j-1][1].as<int>();
                    x2 = ip_result[j][0].as<int>();
                    y2 = ip_result[j][1].as<int>();
                    x_dist = x2 - x1;
                    y_dist = y2 - y1;

                    if (x_dist < x_dim && y_dist < y_dim) {
                        long_centroid = (x1 + x2) / 2;
                        lat_centroid = (y1 + y2) / 2;

                        break;
                    }
                    
                }                     
                long_centroid = x1;
                lat_centroid = y1;     

                int x_off, y_off;
                int rand = std::rand();
                if (x_dim == 1) {
                    x_off = long_centroid;
                } else if (long_centroid < (N - max_x_offset) - 1) {  // point in existing range
                    // possible range [0, x_dim - ((x_dim - longitude)]
                    x_off =rand % std::min((x_dim - (x_dim - long_centroid) + 1), max_x_offset + 1);
                } else {  // out
                    // possible range [longitude - x_dim + 1, min(max_offset, longitude)]
                    int min_off = long_centroid - x_dim + 1;
                    x_off = min_off;
                    if ((std::min(max_x_offset, long_centroid) - min_off) != 0) {
                        x_off += rand % (std::min(max_x_offset, long_centroid) - min_off);
                    }
                }
                //x_off = (int) ScaleToResolutionCeil(ScaleToResolutionCeil((x_off), N, resolution), resolution, N);

                if (y_dim == 1) {
                    y_off = lat_centroid;
                } else if (lat_centroid < (N - max_y_offset) - 1) {
                    y_off = std::rand() % std::min((y_dim - (y_dim - lat_centroid) + 1), max_y_offset + 1);
                } else {  // out
                    int min_off = lat_centroid - y_dim + 1;
                    y_off = min_off;
                    if ((std::min(max_y_offset, lat_centroid) - min_off) != 0) {
                        y_off += std::rand() % (std::min(max_y_offset, lat_centroid) - min_off);
                    }
                }
                //y_off = ScaleToResolutionCeil(ScaleToResolutionCeil((y_off), N, resolution), resolution, N);

                bool item_within_region = false;
                offset_ranges = zerod_ranges;
                for (int j = 0; j < (int) zerod_ranges.size(); j++) {
                    range r = zerod_ranges[j];

                    /*r.x1 = ScaleToResolutionFloor(ScaleToResolutionFloor(r.x1, N, resolution), resolution, N);
                    r.y1 = ScaleToResolutionFloor(ScaleToResolutionFloor(r.y1, N, resolution), resolution, N);
                    r.x2 = ScaleToResolutionCeil(ScaleToResolutionCeil(r.x2, N, resolution), resolution, N) - 1;
                    r.y2 = ScaleToResolutionCeil(ScaleToResolutionCeil(r.y2, N, resolution), resolution, N) - 1;*/

                    r.x1 = r.x1 + x_off; r.x2 += x_off; r.y1 += y_off; r.y2 += y_off;

                    if (r.x1 < 0 || r.x2 >= N || r.y1 < 0 || r.y2 >= N) {
                        //std::cout << "Out of bounds " << r.x1 << ", " << r.y1 << ", " << r.x2 << ", " << r.y2 << " N " << N << std::endl;
                        continue;
                    }
                    if (r.x1 <= long_centroid && long_centroid <= r.x2 && r.y1 <= lat_centroid && lat_centroid <= r.y2) {
                        item_within_region = true;
                    }
                    offset_ranges[j] = r;
                }
                if (item_within_region == false) {
                    //std::cout << "Item not in region" << std::endl;
                    retries++;
                    if (retries > 5 * sample_size) {
                        std::cout << "Quiting GetIPRangeQueries2() after " << 5*sample_size << " tries for a good query set, not in range" << std::endl;
                    }
                    continue;
                }

                long start_ip = 0;
                long end_ip = 0;

                if (range_queries) {
                    // Get all ip addresses in this cell to derive a range, if there are less than 2 ips in this cell, skip
                    // TODO: This does result in sorted ip's, therefore favoring ip's that are close to each other
                    pqxx::work neighbour_ips(C);
                    query = "SELECT ip FROM grouped_ip WHERE (" + ComposeRange(offset_ranges) + ");";
                    neighbour_ips_result = neighbour_ips.exec(query);
                    neighbour_ips.commit();

                    // Query range is accepted, now derive ip range
                    long neighbour_ip = 0;
                    for (auto row : neighbour_ips_result) {
                        neighbour_ip = row[0].as<long>();
                        if (neighbour_ip == item) {
                            continue;
                        }
                        if (item > neighbour_ip) {
                            start_ip = neighbour_ip - std::rand() % (neighbour_ip);
                            end_ip = item + std::rand() % (UINT_MAX - item);
                            break;
                        } else {
                            start_ip = item - std::rand() % (item);
                            end_ip = neighbour_ip + std::rand() % (UINT_MAX - neighbour_ip);
                            break;
                        }
                    }
                    if (start_ip == 0 || end_ip == 0) {
                        // Create random range
                        start_ip = item - std::rand() % (item);
                        end_ip = item + std::rand() % (UINT_MAX - item);
                    }
                } else {
                    start_ip = item;
                    end_ip = item;
                }
                //start_ip = 0; //item - std::rand() % (item);
                //end_ip = UINT_MAX; // item + std::rand() % (UINT_MAX - item);

                // conform ip range to resolution
                //start_ip = ScaleToResolution(ScaleToResolution((long) start_ip, (long) UINT_MAX, resolution), resolution, (long) UINT_MAX);
                //end_ip = ScaleToResolution(ScaleToResolution((long) end_ip, (long) UINT_MAX, resolution), resolution, (long) UINT_MAX);

                // Given new query range with item within range, compute the entire query set data
                offset_vertices = vertices;
                for (size_t j = 0; j < vertices.size(); j++) {
                    offset_vertices[j].first += x_off;
                    offset_vertices[j].second += y_off;
                }

                // Ground truth (needed when one ip is mapped to more cells)
                if (range_queries) {
                    query = "SELECT SUM(count) FROM grouped_ip WHERE (" + ComposeRange(offset_ranges) + ") and (ip >= " + std::to_string(start_ip) + " and ip <= " + std::to_string(end_ip) + ");";
                } else if (timestamp) {
                    query = "SELECT timestamp FROM grouped_ip where ip = " + std::to_string(item) + " LIMIT 1;";
                    pqxx::work wtime(C);
                    pqxx::result timewindow = wtime.exec(query);
                    wtime.commit();
                    if (timewindow.size() == 0 || timewindow[0][0].is_null()) {
                        retries++;
                        if (retries > 10 * sample_size) {
                            std::cout << "Quiting GetIPRangeQueriesFrequfency() after " << 10*sample_size << " tries for a good query set, empty timsetamp" << std::endl;
                            break;
                        }
                        continue;
                    }
                    // Get random number between max possible timestamp and zero
                    int qt = timewindow[0][0].as<int>();
                    std::random_device rd;
                    std::mt19937 gen(rd());
                    std::uniform_int_distribution<> distrib(0, qt);
                    time_window = distrib(gen);
                    query = "SELECT SUM(count) FROM grouped_ip WHERE (" + ComposeRange(offset_ranges) + ") and (ip = " + std::to_string(start_ip) + ") and timestamp >= " + std::to_string(time_window) + ";";
                } else {
                    query = "SELECT SUM(count) FROM grouped_ip WHERE (" + ComposeRange(offset_ranges) + ") and (ip = " + std::to_string(start_ip) + ");";
                }
                pqxx::work w3(C);
                pqxx::result qgroundtruth = w3.exec(query);
                w3.commit();
                if (qgroundtruth.size() == 0 || qgroundtruth[0][0].is_null()) {
                    retries++;
                    if (retries > 10 * sample_size) {
                        std::cout << "Quiting GetIPRangeQueries3() after " << 5*sample_size << " tries for a good query set, empty ground truth" << std::endl;
                        break;
                    }
                    continue;
                }
                int groundtruth = qgroundtruth[0][0].as<int>();

                if (groundtruth < min_query_answer) {
                    std::cout << "Query rejected afer ground truth deemed to low" << std::endl;
                    retries++;
                    if (retries > 10 * sample_size) {
                        std::cout << "Quiting GetIPRangeQueries3() after " << 5*sample_size << " tries for a good query set, low answer" << std::endl;
                    }
                    continue;
                } else {
                    samples++;
                }

                if (timestamp) {
                    query = "SELECT count(ip) FROM iplocations WHERE (" + ComposeRange(offset_ranges) + ") and (timestamp >= " + std::to_string(time_window) + ");";
                } else {
                    query = "SELECT (count) FROM grouped_ip WHERE (" + ComposeRange(offset_ranges) + ");";
                }
                pqxx::work w2(C);
                pqxx::result L1_res = w2.exec(query);
                if (L1_res.size() == 0 || L1_res[0][0].is_null()) { std::cout << "This shouldn't happen" << std::endl; }
                int L1 = L1_res[0][0].as<int>();;
                w2.commit();

                QuerySet q = {offset_ranges, offset_vertices, start_ip, end_ip, groundtruth, L1, dataset_size, x_off, y_off, time_window};
                query_sets.push_back(q);
            }
            C.disconnect();
        }
    } catch (const std::exception &e) {
        std::cerr << "Exception at GetIPInRange: " << e.what() << std::endl;
    }
    return query_sets;
}

std::vector<QuerySet> Postgres::GetIPRangeQueriesCountDistinct(std::vector<range> zerod_ranges, std::vector<std::pair<float, float>> vertices, int N, int max_x_offset, int max_y_offset, int sample_size, int limit, int min_query_answer) {
    // xerod_ranges: 
    // vertices:
    // N: stream size
    // max_x_offset: 
    // max_y_offset:
    // sample_size: number of times same shape is queried on different position.
    // limit:
    // min_query_answer: minimal number of distinct ip's in query range

    
    std::vector<range> offset_ranges;
    std::vector<std::pair<float, float>> offset_vertices;
    std::vector<QuerySet> query_sets = {};
    query_sets.reserve(sample_size);
    try {
        // Todo: Possibly it is not required to open connection each and every time
        pqxx::connection C("dbname = postgres user = postgres password = postgres hostaddr = 127.0.0.1 port = 5432");
        if (C.is_open()) {
            // Get dataset size
            pqxx::work Nquery(C);
            int dataset_size = Nquery.exec("SELECT COUNT(*) FROM iplocations;")[0][0].as<int>();
            Nquery.commit();
            // pqxx::work CDquery(C);
            // int distinct_ips = CDquery.exec("SELECT COUNT(distinct ip) FROM iplocations;")[0][0].as<int>();
            // CDquery.commit();
            // std::cout << "Distinct ips " << distinct_ips << std::endl;

            // Get top N items, their counts and location
            // By creating materialized view first
            pqxx::work mview(C);
            mview.exec("CREATE MATERIALIZED VIEW IF NOT EXISTS grouped_ip_cd AS (SELECT longitude, latitude, count(distinct ip) as dc_ip FROM iplocations GROUP BY longitude, latitude ORDER BY longitude, latitude);");
            mview.exec("CREATE INDEX IF NOT EXISTS mview_btree ON grouped_ip_cd USING btree (longitude, latitude);");
            mview.commit();
            
            // Get all x y locations with at least three distinct ip adresses.
            std::string query = "SELECT longitude, latitude, dc_ip FROM grouped_ip_cd WHERE dc_ip > 3;"; // opt 
            pqxx::work w(C);
            pqxx::result result = w.exec(query);
            w.commit();

            int x_dim = (N - max_x_offset) - 1;
            int y_dim = (N - max_y_offset) - 1;

            int retries = 0;
            int samples = 0;
            int i = -1;
            while (samples < sample_size) {
                i++;
                pqxx::row entry;
                entry = result[(int) i % result.size()];
                long item = entry[0].as<long>();
                long lat = entry[0].as<long>();
                long lon = entry[1].as<long>();

                // Order differently to reduce bias towards some geolocation
                // if (i % 4 == 0) {
                //     query = "SELECT dc_ip FROM grouped_ip_cd WHERE lat = " + std::to_string(lat) + " AND lon = " + std::to_string(lon) + " ORDER BY longitude, latitude DESC;";
                // } else if (i % 4 == 1) {
                //     query = "SELECT dc_ip FROM grouped_ip_cd WHERE lat = " + std::to_string(lat) + " AND lon = " + std::to_string(lon) + " ORDER BY longitude, latitude DESC;";
                    
                //     query = "SELECT longitude, latitude FROM grouped_ip_cd WHERE ip = " + std::to_string(item) + " ORDER BY latitude, longitude DESC;";
                // } else if (i % 4 == 2) {
                //     query = "SELECT longitude, latitude FROM grouped_ip_cd WHERE ip = " + std::to_string(item) + " ORDER BY longitude, latitude ASC;";
                // } else {
                //     query = "SELECT longitude, latitude FROM grouped_ip_cd WHERE ip = " + std::to_string(item) + " ORDER BY latitude, longitude ASC;";
                // }
                pqxx::work ips(C);
                pqxx::result ip_result = ips.exec(query);
                ips.commit();

                int x1, y1, x2, y2, x_dist, y_dist;
                int long_centroid, lat_centroid;
                for (int j = 1; j < (int) ip_result.size(); j++) {
                    x1 = ip_result[j-1][0].as<int>();
                    y1 = ip_result[j-1][1].as<int>();
                    x2 = ip_result[j][0].as<int>();
                    y2 = ip_result[j][1].as<int>();
                    x_dist = x2 - x1;
                    y_dist = y2 - y1;

                    if (x_dist < x_dim && y_dist < y_dim) {
                        long_centroid = (x1 + x2) / 2;
                        lat_centroid = (y1 + y2) / 2;
                        break;
                    }
                }

                int x_off, y_off;
                int rand = std::rand();
                if (x_dim == 1) {
                    x_off = long_centroid;
                } else if (long_centroid < (N - max_x_offset) - 1) {  // point in existing range
                    // possible range [0, x_dim - ((x_dim - longitude)]
                    x_off =rand % std::min((x_dim - (x_dim - long_centroid) + 1), max_x_offset + 1);
                } else {  // out
                    // possible range [longitude - x_dim + 1, min(max_offset, longitude)]
                    int min_off = long_centroid - x_dim + 1;
                    x_off = min_off;
                    if ((std::min(max_x_offset, long_centroid) - min_off) != 0) {
                        x_off += rand % (std::min(max_x_offset, long_centroid) - min_off);
                    }
                }

                if (y_dim == 1) {
                    y_off = lat_centroid;
                } else if (lat_centroid < (N - max_y_offset) - 1) {
                    y_off = std::rand() % std::min((y_dim - (y_dim - lat_centroid) + 1), max_y_offset + 1);
                } else {  // out
                    int min_off = lat_centroid - y_dim + 1;
                    y_off = min_off;
                    if ((std::min(max_y_offset, lat_centroid) - min_off) != 0) {
                        y_off += std::rand() % (std::min(max_y_offset, lat_centroid) - min_off);
                    }
                }

                bool item_within_region = false;
                offset_ranges = zerod_ranges;
                for (int j = 0; j < (int) zerod_ranges.size(); j++) {
                    range r = zerod_ranges[j];
                    r.x1 = r.x1 + x_off; r.x2 += x_off; r.y1 += y_off; r.y2 += y_off;
                    if (r.x1 < 0 || r.x2 >= N || r.y1 < 0 || r.y2 >= N) {
                        std::cout << "Out of bounds " << r.x1 << ", " << r.y1 << ", " << r.x2 << ", " << r.y2 << " N " << N << std::endl;
                        continue;
                    }
                    if (r.x1 <= long_centroid && long_centroid <= r.x2 && r.y1 <= lat_centroid && lat_centroid <= r.y2) {
                        item_within_region = true;
                    }
                    offset_ranges[j] = r;
                }
                if (item_within_region == false) {
                    //std::cout << "Item not in region" << std::endl;
                    retries++;
                    if (retries > 5 * sample_size) {
                        std::cout << "Quiting GetIPRangeQueriesCountDistinct() after " << 5*sample_size << " tries for a good query set, not in range" << std::endl;
                    }
                    continue;
                }

                // Given new query range with item within range, compute the entire query set data
                offset_vertices = vertices;
                for (size_t j = 0; j < vertices.size(); j++) {
                    offset_vertices[j].first += x_off;
                    offset_vertices[j].second += y_off;
                }

                // Ground truth (needed when one ip is mapped to more cells)'
                //query = "SELECT SUM(dc_ip) FROM grouped_ip_cd WHERE (" + ComposeRange(offset_ranges) + ");";

                query = "SELECT Count(distinct ip) FROM iplocations WHERE (" + ComposeRange(offset_ranges) + ");"; //TODO: place back distinct

                pqxx::work w3(C);
                pqxx::result qgroundtruth = w3.exec(query);
                w3.commit();
                if (qgroundtruth.size() == 0 || qgroundtruth[0][0].is_null()) {
                    retries++;
                    if (retries > 5 * sample_size) {
                        std::cout << "Quiting GetIPRangeQueriesCountDistinct() after " << 5*sample_size << " tries for a good query set, empty ground truth" << std::endl;
                    }
                    continue;
                }
                
                int groundtruth = qgroundtruth[0][0].as<int>();
    

                // if (groundtruth < distinct_ips) {
                //     std::cout << "Query found less distinct ips than in dataset, namely: " << groundtruth << "  out of " << distinct_ips<< std::endl;
                // } else {
                //     std::cout << "Query found more distinct ips than in dataset" << std::endl;
                // }
                if (groundtruth < min_query_answer) {
                    std::cout << "Query rejected afer ground truth deemed to low" << std::endl;
                    retries++;
                    if (retries > 5 * sample_size) {
                        std::cout << "Quiting GetIPRangeQueriesCountDistinct() after " << 5*sample_size << " tries for a good query set, low answer" << std::endl;
                    }
                    continue;
                } else {
                    samples++;
                }

                query = "SELECT count(ip) FROM iplocations WHERE (" + ComposeRange(offset_ranges) + ");";
                
                //query = "SELECT SUM(dc_ip) FROM grouped_ip_cd WHERE (" + ComposeRange(offset_ranges) + ");";
                pqxx::work w2(C);
                pqxx::result L1_res = w2.exec(query);
                if (L1_res.size() == 0 || L1_res[0][0].is_null()) { std::cout << "This shouldn't happen" << std::endl; }
                int L1 = L1_res[0][0].as<int>();;
                w2.commit();

                QuerySet q;
                q.ranges = offset_ranges;
                q.vertices = offset_vertices;
                q.offset_x = x_off;
                q.offset_y = y_off;
                q.groundtruth = groundtruth;
                q.N = dataset_size;
                q.L1 = L1;
                query_sets.push_back(q);
            }
            C.disconnect();
        }
    } catch (const std::exception &e) {
        std::cerr << "Exception at GetIPInRange: " << e.what() << std::endl;
    }
    return query_sets;
}


//TODO: This function is not yet finished
std::vector<QuerySet> Postgres::GetIPRangeQueriesMembership(std::vector<range> zerod_ranges, std::vector<std::pair<float, float>> vertices, int N, int max_x_offset, int max_y_offset, int sample_size, int limit, int min_query_answer) {
    // xerod_ranges: 
    // vertices:
    // N: stream size
    // max_x_offset: 
    // max_y_offset:
    // sample_size: number of times same shape is queried on different position.
    // limit:
    // min_query_answer: minimal number of distinct ip's in query range


    std::vector<range> offset_ranges;
    std::vector<std::pair<float, float>> offset_vertices;
    std::vector<QuerySet> query_sets = {};
    query_sets.reserve(sample_size);
    try {
        // Todo: Possibly it is not required to open connection each and every time
        pqxx::connection C("dbname = postgres user = postgres password = postgres hostaddr = 127.0.0.1 port = 5432");
        if (C.is_open()) {
            // Get dataset size
            pqxx::work Nquery(C);
            int dataset_size = Nquery.exec("SELECT COUNT(*) FROM iplocations;")[0][0].as<int>();
            Nquery.commit();
            

            // Get top N items, their counts and location
            // By creating materialized view first
            pqxx::work mview(C);
            mview.exec("CREATE MATERIALIZED VIEW IF NOT EXISTS grouped_ip_mem AS (SELECT longitude, latitude, count(distinct ip) as dc_ip, min(ip) as min_ip, max(ip) as max_ip FROM iplocations GROUP BY longitude, latitude ORDER BY longitude, latitude);");
            mview.exec("CREATE INDEX IF NOT EXISTS mview_btree ON grouped_ip_mem USING btree (longitude, latitude);");
            mview.commit();

            // pqxx::work mview(C);
            // mview.exec("CREATE MATERIALIZED VIEW IF NOT EXISTS distinct_ips AS (SELECT distinct ip as dip FROM iplocations LIMIT 100);");
            // mview.exec("CREATE INDEX IF NOT EXISTS distinct_ips_btree ON distinct_ips USING btree (dip);");
            // mview.commit();
            std::string dist_query = "SELECT distinct ip FROM iplocations LIMIT " + std::to_string(sample_size*5) + ";";
            pqxx::work dist(C);
            pqxx::result distinct_ips = dist.exec(dist_query);
            dist.commit();

            // Get all x y locations with at least three distinct ip adresses.
            // Get x y locations with at least one ip address.
            //std::string query = "SELECT longitude, latitude, dc_ip, min_ip, max_ip FROM grouped_ip_mem GROUP BY longitude, latitude, dc_ip HAVING dc_ip > 3 ORDER BY dc_ip DESC;"; // opt 
            std::string query = "SELECT longitude, latitude, dc_ip, min_ip, max_ip FROM grouped_ip_mem WHERE dc_ip > 3 ORDER BY dc_ip DESC;"; // opt 

            pqxx::work w(C);
            pqxx::result result = w.exec(query);
            w.commit();

            int x_dim = (N - max_x_offset) - 1;
            int y_dim = (N - max_y_offset) - 1;

            int retries = 0;
            int samples = 0;
            int i = -1;
            while (samples < sample_size) {
                i++;
                pqxx::row entry;
                entry = result[(int) i % result.size()];
                //uint item = entry[0].as<uint>();
                long lat = entry[0].as<long>();
                long lon = entry[1].as<long>();
                long min_ip = entry[3].as<long>();
                long max_ip = entry[4].as<long>();
                // Order differently to reduce bias towards some geolocation
                // if (i % 4 == 0) {
                //     query = "SELECT dc_ip FROM grouped_ip_cd WHERE lat = " + std::to_string(lat) + " AND lon = " + std::to_string(lon) + " ORDER BY longitude, latitude DESC;";
                // } else if (i % 4 == 1) {
                //     query = "SELECT dc_ip FROM grouped_ip_cd WHERE lat = " + std::to_string(lat) + " AND lon = " + std::to_string(lon) + " ORDER BY longitude, latitude DESC;";
                    
                //     query = "SELECT longitude, latitude FROM grouped_ip_cd WHERE ip = " + std::to_string(item) + " ORDER BY latitude, longitude DESC;";
                // } else if (i % 4 == 2) {
                //     query = "SELECT longitude, latitude FROM grouped_ip_cd WHERE ip = " + std::to_string(item) + " ORDER BY longitude, latitude ASC;";
                // } else {
                //     query = "SELECT longitude, latitude FROM grouped_ip_cd WHERE ip = " + std::to_string(item) + " ORDER BY latitude, longitude ASC;";
                // }
                pqxx::work ips(C);
                pqxx::result ip_result = ips.exec(query);
                ips.commit();
                
                int x1, y1, x2, y2, x_dist, y_dist;
                int long_centroid, lat_centroid;
                for (int j = 1; j < (int) ip_result.size(); j++) {
                    x1 = ip_result[j-1][0].as<int>();
                    y1 = ip_result[j-1][1].as<int>();
                    x2 = ip_result[j][0].as<int>();
                    y2 = ip_result[j][1].as<int>();
                    x_dist = x2 - x1;
                    y_dist = y2 - y1;

                    if (x_dist < x_dim && y_dist < y_dim) {
                        long_centroid = (x1 + x2) / 2;
                        lat_centroid = (y1 + y2) / 2;
                        break;
                    }
                }

                int x_off, y_off;
                int rand = std::rand();
                if (x_dim == 1) {
                    x_off = long_centroid;
                } else if (long_centroid < (N - max_x_offset) - 1) {  // point in existing range
                    // possible range [0, x_dim - ((x_dim - longitude)]
                    x_off =rand % std::min((x_dim - (x_dim - long_centroid) + 1), max_x_offset + 1);
                } else {  // out
                    // possible range [longitude - x_dim + 1, min(max_offset, longitude)]
                    int min_off = long_centroid - x_dim + 1;
                    x_off = min_off;
                    if ((std::min(max_x_offset, long_centroid) - min_off) != 0) {
                        x_off += rand % (std::min(max_x_offset, long_centroid) - min_off);
                    }
                }

                if (y_dim == 1) {
                    y_off = lat_centroid;
                } else if (lat_centroid < (N - max_y_offset) - 1) {
                    y_off = std::rand() % std::min((y_dim - (y_dim - lat_centroid) + 1), max_y_offset + 1);
                } else {  // out
                    int min_off = lat_centroid - y_dim + 1;
                    y_off = min_off;
                    if ((std::min(max_y_offset, lat_centroid) - min_off) != 0) {
                        y_off += std::rand() % (std::min(max_y_offset, lat_centroid) - min_off);
                    }
                }

                bool item_within_region = false;
                offset_ranges = zerod_ranges;
                for (int j = 0; j < (int) zerod_ranges.size(); j++) {
                    range r = zerod_ranges[j];
                    r.x1 = r.x1 + x_off; r.x2 += x_off; r.y1 += y_off; r.y2 += y_off;
                    if (r.x1 < 0 || r.x2 >= N || r.y1 < 0 || r.y2 >= N) {
                        std::cout << "Out of bounds " << r.x1 << ", " << r.y1 << ", " << r.x2 << ", " << r.y2 << " N " << N << std::endl;
                        continue;
                    }
                    if (r.x1 <= long_centroid && long_centroid <= r.x2 && r.y1 <= lat_centroid && lat_centroid <= r.y2) {
                        item_within_region = true;
                    }
                    offset_ranges[j] = r;
                }
                if (item_within_region == false) {
                    //std::cout << "Item not in region" << std::endl;
                    retries++;
                    if (retries > 5 * sample_size) {
                        std::cout << "Quiting GetIPRangeQueriesCountDistinct() after " << 5*sample_size << " tries for a good query set, not in range" << std::endl;
                    }
                    continue;
                }

                // Given new query range with item within range, compute the entire query set data
                offset_vertices = vertices;
                for (size_t j = 0; j < vertices.size(); j++) {
                    offset_vertices[j].first += x_off;
                    offset_vertices[j].second += y_off;
                }

                // Ground truth (needed when one ip is mapped to more cells)'
                
                pqxx::row dist_entry = distinct_ips[(int) i % distinct_ips.size()];
                long queried_ip_long = dist_entry[0].as<long>();
                
                std::string queried_ip = std::to_string(queried_ip_long);
                query = "SELECT EXISTS(SELECT 1 FROM iplocations WHERE (" + ComposeRange(offset_ranges) + ") AND (ip = " + queried_ip + "));";
                pqxx::work w3(C);
                pqxx::result qgroundtruth = w3.exec(query);
                w3.commit();
    
                std::string s = qgroundtruth[0][0].as<std::string>();
                int groundtruth;
                if (s == "t") {
                    groundtruth = 1;
                } else {
                    groundtruth = 0;
                }
                //int groundtruth;// = result[0][0].as<int>();
                //std::cout << "ground truth is of form " <<qgroundtruth.asString() << std::endl;
                // if (groundtruth < min_query_answer) {
                //     std::cout << "Query rejected afer ground truth deemed to low" << std::endl;
                //     retries++;
                //     if (retries > 5 * sample_size) {
                //         std::cout << "Quiting GetIPRangeQueriesCountDistinct() after " << 5*sample_size << " tries for a good query set, low answer" << std::endl;
                //     }
                //     continue;
                // } else {
                samples++;
                // }
                //query = "SELECT EXISTS(SELECT 1 FROM iplocations WHERE (" + ComposeRange(offset_ranges) + ") AND (ip = " + queried_ip + "));";
                // pqxx::work w2(C);
                // pqxx::result L1_res = w2.exec(query);
                // if (L1_res.size() == 0 || L1_res[0][0].is_null()) { std::cout << "This shouldn't happen" << std::endl; }
                // int L1 = L1_res[0][0].as<int>();;
                // w2.commit();
                //TODO: CHECK queryset content and what item is.
                QuerySet q;
                q.ranges = offset_ranges;
                q.vertices = offset_vertices;
                q.offset_x = x_off;
                q.offset_y = y_off;
                q.item = queried_ip_long;
                q.groundtruth = groundtruth;
                q.N = dataset_size;
                query_sets.push_back(q);
            }
            C.disconnect();
        }
    } catch (const std::exception &e) {
        std::cerr << "Exception at GetIPInRange: " << e.what() << std::endl;
    }
    return query_sets;
}

std::vector<QuerySet> Postgres::GetIPRangeQueriesL2(std::vector<range> zerod_ranges, std::vector<std::pair<float, float>> vertices, int N, int max_x_offset, int max_y_offset, int sample_size, int min_query_answer) {
    std::vector<range> offset_ranges;
    std::vector<std::pair<float, float>> offset_vertices;
    std::vector<QuerySet> query_sets = {};
    query_sets.reserve(sample_size);
    pqxx::result neighbour_ips_result;
    try {
        // Todo: Possibly it is not required to open connection each and every time
        pqxx::connection C("dbname = postgres user = postgres password = postgres hostaddr = 127.0.0.1 port = 5432");
        if (C.is_open()) {
            // Get dataset size
            pqxx::work Nquery(C);
            int dataset_size = Nquery.exec("SELECT COUNT(*) FROM iplocations;")[0][0].as<int>();
            Nquery.commit();

            // Get top N items, their counts and location
            // By creating materialized view first
            pqxx::work mview(C);
            mview.exec("CREATE MATERIALIZED VIEW IF NOT EXISTS grouped_ip AS (SELECT ip, longitude, latitude, count(ip) FROM iplocations GROUP BY ip, longitude, latitude ORDER BY ip, longitude, latitude);");
            mview.exec("CREATE INDEX IF NOT EXISTS mview_btree ON grouped_ip USING btree (ip, longitude, latitude);");
            mview.commit();
            
            // Get all ip address that occur at atleast two locations
            std::string query = "SELECT ip, count(ip) FROM grouped_ip GROUP BY ip HAVING count(ip) >= 2 ORDER BY count(ip) DESC LIMIT " + std::to_string(sample_size * 5) + ";"; // opt 
            pqxx::work w(C);
            pqxx::result result = w.exec(query);
            w.commit();

            int x_dim = (N - max_x_offset) - 1;
            int y_dim = (N - max_y_offset) - 1;
            int cmax_x_offset = max_x_offset; //ScaleToResolutionFloor(ScaleToResolutionFloor(N - x_dim, N, resolution), resolution, N);
            int cmax_y_offset = max_y_offset; //ScaleToResolutionFloor(ScaleToResolutionFloor(N - y_dim, N, resolution), resolution, N);

            int retries = 0;
            int samples = 0;
            int i = -1;
            while (samples < sample_size) {
                i++;
                pqxx::row entry;
                entry = result[(int) std::rand() /*i*/ % result.size()];
                
                long item = entry[0].as<long>();

                // Order differently to reduce bias towards some geolocation
                if (i % 4 == 0) {
                    query = "SELECT longitude, latitude FROM grouped_ip WHERE ip = " + std::to_string(item) + " ORDER BY longitude, latitude DESC;";
                } else if (i % 4 == 1) {
                    query = "SELECT longitude, latitude FROM grouped_ip WHERE ip = " + std::to_string(item) + " ORDER BY latitude, longitude DESC;";
                } else if (i % 4 == 2) {
                    query = "SELECT longitude, latitude FROM grouped_ip WHERE ip = " + std::to_string(item) + " ORDER BY longitude, latitude ASC;";
                } else {
                    query = "SELECT longitude, latitude FROM grouped_ip WHERE ip = " + std::to_string(item) + " ORDER BY latitude, longitude ASC;";
                }
                pqxx::work ips(C);
                pqxx::result ip_result = ips.exec(query);
                ips.commit();

                int x1, y1, x2, y2, x_dist, y_dist;
                int long_centroid, lat_centroid;
                for (int j = 1; j < (int) ip_result.size(); j++) {
                    x1 = ip_result[j-1][0].as<int>();
                    y1 = ip_result[j-1][1].as<int>();
                    x2 = ip_result[j][0].as<int>();
                    y2 = ip_result[j][1].as<int>();
                    x_dist = x2 - x1;
                    y_dist = y2 - y1;

                    if (x_dist < x_dim && y_dist < y_dim) {
                        long_centroid = (x1 + x2) / 2;
                        lat_centroid = (y1 + y2) / 2;

                        break;
                    }
                    
                }                     
                long_centroid = x1;
                lat_centroid = y1;     

                int x_off, y_off;
                int rand = std::rand();
                if (x_dim == 1) {
                    x_off = long_centroid;
                } else if (long_centroid < (N - max_x_offset) - 1) {  // point in existing range
                    // possible range [0, x_dim - ((x_dim - longitude)]
                    x_off =rand % std::min((x_dim - (x_dim - long_centroid) + 1), max_x_offset + 1);
                } else {  // out
                    // possible range [longitude - x_dim + 1, min(max_offset, longitude)]
                    int min_off = long_centroid - x_dim + 1;
                    x_off = min_off;
                    if ((std::min(max_x_offset, long_centroid) - min_off) != 0) {
                        x_off += rand % (std::min(max_x_offset, long_centroid) - min_off);
                    }
                }
                //x_off = (int) ScaleToResolutionCeil(ScaleToResolutionCeil((x_off), N, resolution), resolution, N);

                if (y_dim == 1) {
                    y_off = lat_centroid;
                } else if (lat_centroid < (N - max_y_offset) - 1) {
                    y_off = std::rand() % std::min((y_dim - (y_dim - lat_centroid) + 1), max_y_offset + 1);
                } else {  // out
                    int min_off = lat_centroid - y_dim + 1;
                    y_off = min_off;
                    if ((std::min(max_y_offset, lat_centroid) - min_off) != 0) {
                        y_off += std::rand() % (std::min(max_y_offset, lat_centroid) - min_off);
                    }
                }
                //y_off = ScaleToResolutionCeil(ScaleToResolutionCeil((y_off), N, resolution), resolution, N);

                bool item_within_region = false;
                offset_ranges = zerod_ranges;
                for (int j = 0; j < (int) zerod_ranges.size(); j++) {
                    range r = zerod_ranges[j];

                    /*r.x1 = ScaleToResolutionFloor(ScaleToResolutionFloor(r.x1, N, resolution), resolution, N);
                    r.y1 = ScaleToResolutionFloor(ScaleToResolutionFloor(r.y1, N, resolution), resolution, N);
                    r.x2 = ScaleToResolutionCeil(ScaleToResolutionCeil(r.x2, N, resolution), resolution, N) - 1;
                    r.y2 = ScaleToResolutionCeil(ScaleToResolutionCeil(r.y2, N, resolution), resolution, N) - 1;*/

                    r.x1 = r.x1 + x_off; r.x2 += x_off; r.y1 += y_off; r.y2 += y_off;

                    if (r.x1 < 0 || r.x2 >= N || r.y1 < 0 || r.y2 >= N) {
                        //std::cout << "Out of bounds " << r.x1 << ", " << r.y1 << ", " << r.x2 << ", " << r.y2 << " N " << N << std::endl;
                        continue;
                    }
                    if (r.x1 <= long_centroid && long_centroid <= r.x2 && r.y1 <= lat_centroid && lat_centroid <= r.y2) {
                        item_within_region = true;
                    }
                    offset_ranges[j] = r;
                }
                if (item_within_region == false) {
                    //std::cout << "Item not in region" << std::endl;
                    retries++;
                    if (retries > 5 * sample_size) {
                        std::cout << "Quiting GetIPRangeQueries2() after " << 5*sample_size << " tries for a good query set, not in range" << std::endl;
                    }
                    continue;
                }

                // Given new query range with item within range, compute the entire query set data
                offset_vertices = vertices;
                for (size_t j = 0; j < vertices.size(); j++) {
                    offset_vertices[j].first += x_off;
                    offset_vertices[j].second += y_off;
                }

                // Ground truth (needed when one ip is mapped to more cells)
                //query = "SELECT SUM(count * count) FROM (SELECT ip, count FROM grouped_ip WHERE (" + ComposeRange(offset_ranges) + ") GROUP BY ip, count);";
                query = "SELECT SUM(count * count) FROM (SELECT ip, count(*) FROM iplocations WHERE (" + ComposeRange(offset_ranges) + ") GROUP BY ip);";
                //SELECT SUM(cnt*cnt) FROM (SELECT ip, SUM(val) as cnt FROM tmp2 WHERE x,y in range GROUP BY IP);
                pqxx::work w3(C);
                pqxx::result qgroundtruth = w3.exec(query);
                w3.commit();
                if (qgroundtruth.size() == 0 || qgroundtruth[0][0].is_null()) {
                    retries++;
                    if (retries > 10 * sample_size) {
                        std::cout << "Quiting GetIPRangeQueries3() after " << 5*sample_size << " tries for a good query set, empty ground truth" << std::endl;
                        break;
                    }
                    continue;
                }
                long groundtruth = qgroundtruth[0][0].as<long>();

                if (groundtruth < min_query_answer) {
                    std::cout << "Query rejected afer ground truth deemed to low" << std::endl;
                    retries++;
                    if (retries > 10 * sample_size) {
                        std::cout << "Quiting GetIPRangeQueries3() after " << 5*sample_size << " tries for a good query set, low answer" << std::endl;
                    }
                    continue;
                } else {
                    samples++;
                }

                query = "SELECT SUM(count) FROM grouped_ip WHERE (" + ComposeRange(offset_ranges) + ");";
                pqxx::work w2(C);
                pqxx::result L1_res = w2.exec(query);
                if (L1_res.size() == 0 || L1_res[0][0].is_null()) { std::cout << "This shouldn't happen" << std::endl; }
                int L1 = L1_res[0][0].as<int>();;
                w2.commit();

                QuerySet q = {offset_ranges, offset_vertices, 0, 0, groundtruth, L1, dataset_size, x_off, y_off, 0};
                query_sets.push_back(q);
            }
            C.disconnect();
        }
    } catch (const std::exception &e) {
        std::cerr << "Exception at GetIPInRange: " << e.what() << std::endl;
    }
    return query_sets;
}
