//
// Created by yangyueji on 8/29/19.
//

#include "RandSim.h"
#include <queue>
#include <iostream>
#include <random>
using namespace std;
void RandSim::randomWalkGeoSampling(const GraphManager &gm,
                                    std::unordered_map<uint32_t, double> &source2weight,
                                    std::unordered_map<uint32_t, double> &node2score,
                                    double teleprob) {

    int sampleNum = 1000000;
    vector<uint32_t> sources;
    vector<double> initProb;
    for (auto kv_pair : source2weight) {
        sources.push_back(kv_pair.first);
        initProb.push_back(kv_pair.second);
    }
    default_random_engine generator;
    discrete_distribution<int> d(initProb.begin(), initProb.end());
    std::uniform_int_distribution<int> unid(1,100);
    //srand(time(NULL));
    for (int i = 0; i < sampleNum; ++i) {
        int s = sources[d(generator)];
        //generate a random walk
        int r = unid(generator);
        while (r > 100 * teleprob) {
            //stop
            int edge_start = gm.nodes[s];
            int adj_sz = gm.nodes[s + 1] - edge_start;
            uniform_int_distribution<int> tmp_uni(0, adj_sz - 1);
            int pos = tmp_uni(generator);
            s = extractLow32bits(gm.edges[edge_start + pos]);
            r = unid(generator);
        }
        if (node2score.count(s)) {
            node2score[s]++;
        } else {
            node2score[s] = 1;
        }
    }

    //vector<pair<double, uint32_t>> scoreNidPair;
    for (auto & kv_pair : node2score) {
        kv_pair.second /= sampleNum;
     //   scoreNidPair.emplace_back(make_pair(kv_pair.second, kv_pair.first));
    }
//
//    sort(scoreNidPair.begin(), scoreNidPair.end());
//    reverse(scoreNidPair.begin(), scoreNidPair.end());
//    int printSz = 5;
//    for (int k = 0; k < printSz && k < scoreNidPair.size(); ++k) {
//        string iri = gm.nid2IRI[scoreNidPair[k].second];
//        printf("IRI = %s, score = %lf\n", iri.c_str(), scoreNidPair[k].first);
//    }


//    exit(11);
    // print BFS neighborhood
//    for (auto s : sources) {
////        cout << "--------------------------------------------------\n";
//        queue<uint32_t> q;
//        q.push(s);
//
//        unordered_set<uint32_t> visited;
//        for (int i = 0; i < printSz && !q.empty(); ++i) {
//            int cnode = q.front();
//            q.pop();
//            if (visited.count(cnode)) continue;
//            visited.insert(cnode);
//            string iri = gm.nid2IRI[cnode];
//            double score = 0;
//            if (node2score.count(cnode))
//                score = node2score[cnode];
////            printf("IRI = %s, score = %lf\n", iri.c_str(), score);
//
//            //expansion
//            int edge_start = gm.nodes[cnode];
//            int adj_sz = gm.nodes[cnode + 1] - edge_start;
//            for (int j = 0; j < adj_sz; ++j) {
//                auto end_node = extractLow32bits(gm.edges[edge_start + j]);
//                if (visited.count(end_node)) continue;
//                q.push(end_node);
//            }
//        }
//    }
}