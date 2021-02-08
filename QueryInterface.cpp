//
// Created by yangyueji on 10/8/19.
//

//#include <fstream>
#include <iostream>
#include <dirent.h>
#include "QueryInterface.h"
#include "FMiner.h"

using namespace std;


void QueryInterface::getQueryIRIPairs(GraphManager &gm, const char *filename, std::string nidQueryFile) {
    ifstream ifs(filename);
    if (!ifs.is_open()) {
        cout << "Query Node IRI File not exist!" << endl;
        exit(123);
    }
    queryNidPairs.clear();
    queryIRIPairs.clear();
    unordered_set<string> IRIs;
    unordered_map<string, uint32_t> iri2Nid;
    string line;
    while (getline(ifs, line)) {
        if (line.empty()) continue;
        istringstream iss(line);
        string target, context;
        iss >> target >> context;
        IRIs.insert(target);
        IRIs.insert(context);
        queryIRIPairs.emplace_back(make_pair(target, context));
    }
    ifs.close();
    cout << "Loading query IRI finish. # queries = " << queryIRIPairs.size() << endl;
    string prefix = "https://www.wikidata.org/wiki/";
    for (uint32_t i = 0; i < gm.node_sz; ++i) {
        if (iri2Nid.size() == IRIs.size()) break;
        string subs = gm.nid2IRI[i].substr(prefix.size() + 1);
        if (IRIs.count(subs)) {
            iri2Nid[subs] = i;
        }
    }

    for (auto &p : queryIRIPairs) {
        queryNidPairs.emplace_back(make_pair(iri2Nid[p.first], iri2Nid[p.second]));
    }
    cout << "Find Node IRI w.r.t. IRIs. # queries = " << queryIRIPairs.size() << endl;


    ofstream ofs(nidQueryFile);
    for (auto &p : queryNidPairs) {
        ofs << to_string(p.first) + " " + to_string(p.second) + "\n";
    }
    ofs.close();
}


void QueryInterface::getQueryNidPairs(GraphManager &gm, const char *filename) {
    ifstream ifs(filename);
    if (!ifs.is_open()) {
        cout << "Query Node Ids File not exist!" << endl;
        exit(123);
    }
    queryNidPairs.clear();
    unordered_set<string> IRIs;
    string line;
    while (getline(ifs, line)) {
        if (line.empty() || line[0] == '#') continue;
        istringstream iss(line);
        string target, context;
        iss >> target >> context;
        //if (target == context) continue;
        IRIs.insert(target);
        IRIs.insert(context);


        if (target[0] == 'Q' && context[0] == 'Q') {
            if (!gm.iri2Nid.count(target) || !gm.iri2Nid.count(context)) continue;
            queryNidPairs.emplace_back(make_pair(gm.iri2Nid[target], gm.iri2Nid[context]));
        } else
            queryNidPairs.emplace_back(make_pair((uint32_t) stoi(target), (uint32_t) stoi(context)));
    }
    ifs.close();
    //cout << "Loading query node ids finish." << endl;
}

void QueryInterface::processAllQueries(GraphManager &gm, std::string nidQueryFile,
                                       int _top_k_pattern,
                                       int _length_limit,
                                       bool _reverse_search,
                                       bool _optDFS) {
    //get queries
    getQueryNidPairs(gm, nidQueryFile.c_str());

    if (queryNidPairs.empty()) {
        cout << "No query loaded!" << endl;
        return;
    }

    //loop through queries for processing
    toBePreserved.clear();
    double alpha = 0.1;
    double avg_time = 0, effecCount = 0;
    for (int i = 0; i < queryNidPairs.size(); ++i) {

        uint32_t entityInterest = queryNidPairs[i].first;
        uint32_t contextNid = queryNidPairs[i].second;

        std::unordered_map<uint32_t, double> source2weight;
        std::unordered_map<uint32_t, double> node2score;
        RandSim rs;
        source2weight[contextNid] = 1;
        rs.randomWalkGeoSampling(gm, source2weight, node2score, 0.85);

        ExceptionalFactMining efm;
        efm.path_length_limit = _length_limit; // default 4
        efm.num_peerEntity = 5; //significant value
        efm.global_topk = 20;
        efm.top_pattern_k = _top_k_pattern; // used for top-k search, default 20
        efm.degree_threshold = 10000;
//        efm.degree_threshold = 1000;
        efm.nodeScore = node2score;

        efm.contextNodes.insert(contextNid);

        //efm.contextNodes.insert(gm.iri2Nid["Q6581097"]);
        printf("Node ID pair: (%d, %d), degree: (%d, %d)\n", entityInterest, contextNid,
               gm.nodes[entityInterest + 1] - gm.nodes[entityInterest],
               gm.nodes[contextNid + 1] - gm.nodes[contextNid]);
        printf("Entity Interest: %s, score = %.8lf\n", gm.nid2IRI[entityInterest].c_str(),
               efm.nodeScore[entityInterest]);

        if (!node2score.count(entityInterest)) continue;

        efm.relevancy_limit = efm.nodeScore[entityInterest] * alpha;


//        efm.relevancy_limit = 0;
//        efm.oneHopNeighborFind(gm, entityInterest);

        efm.relevancy_limit = 0;
        efm.global_topk = 100;
        efm.enumPathPatternBlind(gm, entityInterest);
//        printf("#!#$#InnerTime=%.6lf\n", efm.InnerProcessingTime);

//        efm.naiveDFS = !_optDFS;
//
//        if (_reverse_search)
//            efm.enumPathPatternWithReverseExpansion_DFS_Reverse(gm, entityInterest, h_num);
//        else
//            efm.enumPathPatternWithReverseExpansion_naive(gm, entityInterest, h_num);

        double gap = 0, s_score = 0, e_score = 0;
        auto num_collected = efm.topk_relevant_pq.size();
        if (!efm.topk_relevant_pq.empty()) {
            printf("The scores are: ");
            s_score = efm.topk_relevant_pq.top().patternRelevance;
            while (!efm.topk_relevant_pq.empty()) {
                e_score = efm.topk_relevant_pq.top().patternRelevance;
                printf("%.8lf, ", e_score);
                efm.topk_relevant_pq.pop();
            }
            gap = e_score - s_score;
        }
        printf("Before top-k PPEvaluation time: %.3lf\n", efm.backwardTimeBeforetopkCollected);
        printf("Dynamic Index Time,  Pattern evaluation time, Traversal time, Total inner time, number collected:%ld, Gap:%.8lf\n#!#%.3lf,%.3lf,%.3lf,%.3lf\n",
               num_collected, gap,
               efm.reverseExpansionTime, efm.backwardFromNodeTime,
               efm.InnerProcessingTime - efm.backwardFromNodeTime - efm.reverseExpansionTime,
               efm.InnerProcessingTime);

//        cout << "--------------------------------" << endl;
        if (!efm.res_pq.empty()) {
            toBePreserved.push_back(queryNidPairs[i]);
            avg_time += efm.InnerProcessingTime;
            effecCount++;
        }
        fflush(stdout);
    }

    printf("Total avg time = %.6lf ms. Effect Count = %lf\n", avg_time / effecCount, effecCount);
}





void QueryInterface::processAllFiles(GraphManager &gm, const char *dir_name, const char *dir2write) {
    DIR *dir;
    struct dirent *ent;
    int count = 0;
    if ((dir = opendir(dir_name)) != NULL) {
        /* print all the files and directories within directory */
        while ((ent = readdir(dir)) != NULL) {
//            printf ("%s\n", ent->d_name);
            string filePath(dir_name);
            string fname(ent->d_name);
            if (fname.find(".txt") == string::npos) continue;
            if (count % 500 == 0)
                cout << " ------------------- " << count << endl;
            count++;
            filePath.append("/" + fname);
            processAllQueries(gm, filePath);
            if (toBePreserved.empty()) continue;
            string outputDir(dir2write);
            outputDir.append("/" + fname);
            write2file(outputDir);
        }
        closedir(dir);
    } else {
        /* could not open directory */
        perror("Cannot open the dir");
        exit(EXIT_FAILURE);
    }

}

std::string get_env1( const std::string & var ) {
    const char * val = std::getenv( var.c_str() );
    if ( val == nullptr ) { // invalid to assign nullptr to std::string
        return "";
    }
    else {
        return val;
    }
}

void QueryInterface::process_all_queries_fminer(GraphManager &gm, std::string nidQueryFile,
                                                int _top_k_pattern,
                                                int _length_limit,
                                                bool delay_high_degree,
                                                bool use_memo,
                                                bool size_bound_opt_on) {
    printf("\n\n-----New Line-----\n\n");
    string INFO = get_env1("INFO");

    //get queries
    getQueryNidPairs(gm, nidQueryFile.c_str());

    if (queryNidPairs.empty()) {
        cout << "No query loaded!" << endl;
        return;
    }

    //loop through queries for processing
    toBePreserved.clear();
    double avg_time = 0, effective_count = 0;
    int test_development_num = 1000; // for dev purpose
    if (queryNidPairs.size() < test_development_num) {
        test_development_num = queryNidPairs.size();
    }

    double ttl_BFS_time = 0;
    double ttl_pattern_eval_time = 0;
    double ttl_node_type_time = 0;
    double ttl_cof_time =0;

    int larger100s = 0;

    for (int i = 0; i < queryNidPairs.size() && i < test_development_num; ++i) {

        uint32_t entityInterest = queryNidPairs[i].first;
        uint32_t contextNid = queryNidPairs[i].second;

        std::unordered_map<uint32_t, double> source2weight;
        std::unordered_map<uint32_t, double> node2score;
        RandSim rs;
        source2weight[contextNid] = 1;
        rs.randomWalkGeoSampling(gm, source2weight, node2score, 0.85);

        FMiner fmr;

#ifdef ServerRun
        if (INFO == "1")
            fmr.verbose = true;
        fmr.path_length_limit = 4; // default 4
        fmr.num_peer_entities = 20; //significant value
        fmr.pattern_topk_num = 20; // top-k patterns
        fmr.cof_topk_num = 20; // top-k cofs per pattern

        fmr.degree_threshold = 10000;

        // Configure the optimization options
        fmr.delay_highdegree_expansion = delay_high_degree;
        fmr.use_memo = use_memo;
        fmr.size_bound_opt_on = size_bound_opt_on;

#else
        fmr.verbose = true;
        fmr.path_length_limit = _length_limit; // default 4
        fmr.num_peer_entities = 2; //significant value
        fmr.pattern_topk_num = 10; // top-k patterns
        fmr.cof_topk_num = 5; // top-k cofs per pattern

        fmr.degree_importance_limit = 2;
        fmr.attribute_occurrence_ratio_limit = 0.1;

        fmr.degree_threshold = 10000;
#endif

        fmr.node_score = node2score;
        fmr.context_nodes.insert(contextNid);

//        printf("Node ID pair: (%d, %d), degree: (%d, %d)\t", entityInterest, contextNid,
//               gm.nodes[entityInterest + 1] - gm.nodes[entityInterest],
//               gm.nodes[contextNid + 1] - gm.nodes[contextNid]);
//        printf("Entity Interest: %s, score = %.8lf\n", gm.nid2IRI[entityInterest].c_str(),
//               fmr.node_score[entityInterest]);

        if (!node2score.count(entityInterest)) continue;

        fmr.FMiner_biBFS(gm, entityInterest);
//        printf("###%d,%.3lf,%.3lf,%.3lf\n", i, fmr.BFS_expansion_time, fmr.pattern_evaluation_time, fmr.COF_extraction_time);
        if (i % 10 == 0)
            fflush(stdout);
//        printf("Time Profiling: BFS_time=%.3lf ms, Pattern_eval_time=%.3lf ms (%d evaluated, %d inserted, %d uninserted, "
//               "(path_score_time=%.3lf ms, node_type_time=%.3lf ms, pattern_enum_time=%.3lf ms, matching_node_set_time=%.3lf ms)), COF_extract_time=%.3lf ms.\n",
//               fmr.BFS_expansion_time,
//               fmr.pattern_evaluation_time, fmr.pattern_evaluation_num, fmr.inserted_pattern_num, fmr.uninserted_pattern_num,
//               fmr.pattern_eval_path_score_time, fmr.pattern_eval_node_type_time, fmr.pattern_enumeration_time,
//               fmr.matching_node_set_call_time,
//               fmr.COF_extraction_time);

        if (fmr.pattern_evaluation_time > 100000) {
            larger100s++;
        }

        ttl_BFS_time += fmr.BFS_expansion_time;
        ttl_pattern_eval_time += fmr.pattern_evaluation_time;
        ttl_node_type_time += fmr.pattern_eval_node_type_time;
        ttl_cof_time += fmr.COF_extraction_time;
        effective_count++;
    }

    printf("AVG:\nBFS_time=%.3lf ms, pattern_eval=%.3lf ms, node_type_time=%.3lf ms, cof_time=%.3lf ms\n",
           ttl_BFS_time/effective_count, ttl_pattern_eval_time/effective_count,
           ttl_node_type_time/effective_count,
           ttl_cof_time/effective_count);
    printf("Effective count= %d, larger than 100s = %d\n", (int) effective_count, larger100s);

    printf("\n\n-----New Line-----\n\n");
    fflush(stdout);
}