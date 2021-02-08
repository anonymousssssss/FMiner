//
// Created by yangyueji on 10/8/19.
//

#ifndef PATHPATTERNMINING_QUERYINTERFACE_H
#define PATHPATTERNMINING_QUERYINTERFACE_H

#include "GraphManager.h"
#include "ExceptionalFactMining.h"
#include "RandSim.h"
#include "settings.h"
#include <fstream>




class QueryInterface {
public:
    int h_num = 2;

    std::vector<std::pair<std::string, std::string>> queryIRIPairs;
    std::vector<std::pair<uint32_t, uint32_t>> queryNidPairs;
    std::vector<std::pair<uint32_t, uint32_t>> toBePreserved;

    ///Assume the file contains a pair of entity IRIs in each line, separated by space.
    void getQueryIRIPairs(GraphManager &gm, const char *filename,
                          std::string nidQueryFile = "../../testData/querySinglePair.txt");

//    std::string nidQueryFile = "../../testData/queryNodeIdPairs.txt";

    //test single pair: Q9616, Q253414. The respective node id is 6737 67205


    void getQueryNidPairs(GraphManager &gm, const char *filename);


    /// In-Loop processing for all queries.
    void processAllQueries(GraphManager &gm, std::string nidQueryFile = "../../testData/querySinglePair.txt",
                           int top_k_pattern = 20,
                           int length_limit = 4,
                           bool reverse_search = true,
                           bool optDFS = true);

    /// Used for new fminer algorithm.
    void process_all_queries_fminer(GraphManager &gm, std::string nidQueryFile = "../../testData/querySinglePair.txt",
                                    int top_k_pattern = 20,
                                    int length_limit = 4,
                                    bool delay_high_degree = true,
                                    bool use_memo = true,
                                    bool size_bound_opt_on = true);

    void processAllFiles(GraphManager &gm, const char *dir_name, const char *path2write);


    inline void write2file(const std::string &fpath) {
        std::ofstream ofs(fpath);
        for (auto &p : toBePreserved) {
            ofs << p.first << " " << p.second << "\n";
        }

        ofs.clear();
    }

};


#endif //PATHPATTERNMINING_QUERYINTERFACE_H
