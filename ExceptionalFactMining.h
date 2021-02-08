//
// Created by yangyueji on 9/13/19.
//

#ifndef PATHPATTERNMINING_EXCEPTIONALFACTMINING_H
#define PATHPATTERNMINING_EXCEPTIONALFACTMINING_H

#include <vector>
#include "GraphManager.h"
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include <queue>
#include <iostream>

typedef std::pair<double, std::string> scoreFact;
#define TimeLimit 300000

struct result_struct {
    double excepScore = 0;
    uint32_t fact_edge;
    uint32_t fact_nodeVal;
    uint32_t current_node;
    std::string pattern_sig;
    double patternRelevance;

    uint32_t mostFrequentValue;
    double freqMFV;
    int current_count = 0;
    std::unordered_set<uint32_t> peer_entities;
    std::vector<std::string> nodeSeqLabels;

    result_struct(double es, uint32_t fe, uint32_t fn, uint32_t cn, std::string &ps, double pr,
                  std::unordered_set<uint32_t> &peers) {
        excepScore = es;
        fact_edge = fe;
        fact_nodeVal = fn;
        current_node = cn;
        pattern_sig = ps;
        patternRelevance = pr;
        peer_entities = peers;
    }

    result_struct(double es, uint32_t fe, uint32_t fn, uint32_t cn, std::string &ps, double pr) {
        excepScore = es;
        fact_edge = fe;
        fact_nodeVal = fn;
        current_node = cn;
        pattern_sig = ps;
        patternRelevance = pr;
    }

    inline void set_node_labels(GraphManager &gm, uint32_t target, std::vector<std::string> &reversedLabels) {
        std::string _l;
        gm.getOneLabelRandom(target, _l);
        nodeSeqLabels.push_back(_l);
        for (int i = 0; i < reversedLabels.size(); ++i) {
            nodeSeqLabels.push_back(reversedLabels[i]);
        }
    }


};

class res_comp {
public:
    bool operator()(const result_struct &r1, const result_struct &r2) {
        return r1.excepScore > r2.excepScore;
    }
};

class res_relevance_comp {
public:
    bool operator()(const result_struct &r1, const result_struct &r2) {
        return r1.patternRelevance > r2.patternRelevance;
    }
};




typedef std::priority_queue<result_struct, std::vector<result_struct>, res_comp> Res_PQ;
typedef std::priority_queue<result_struct, std::vector<result_struct>, res_relevance_comp> Relevant_PQ;



class ExceptionalFactMining {
public:
    int num_peerEntity = 10;
    int path_length_limit = 4;

    double relevancy_limit = 0.000001;
    double degree_threshold = 10000;

    double backwardFromNodeTime = 0;
    double emitResTime = 0;
    double reverseExpansionTime = 0;

    double backwardTimeBeforetopkCollected = 0;

    double InnerProcessingTime = 0;

//    double node_proximity_threshold = 0.000001;
    ///For ease of access
    std::unordered_map<uint32_t, double> nodeScore;
    std::unordered_set<uint32_t> contextNodes;
    std::unordered_set<uint32_t> highDegreeContextNodes;

    int opt_1 = 0;
    int opt_2 = 0;

    double base_e = exp(1);
    //double base_e = 1;
    int global_topk = 10;

    int top_pattern_k = 20;
    Res_PQ res_pq;

    // top-k relevant pattern priority queue
    Relevant_PQ topk_relevant_pq;

    Res_PQ diverseTopk;

    std::unordered_map<uint32_t, double> nid2largestNeighborScore;
//    std::unordered_map<uint32_t, double> eid2largestScore;
    std::vector<std::unordered_map<uint32_t, std::unordered_map<std::string, double>>> hopNid2Parent;
    std::vector<std::unordered_map<uint32_t, std::vector<std::pair<double, std::string>>>> hopNid2Parent_vec;

//    std::unordered_map<std::string, std::unordered_set<uint32_t>> pattern2leftnode;

//    std::unordered_map<std::string, std::vector<std::unordered_set<uint32_t>>> pattern2rightFrontiers;
    std::unordered_map<std::string, std::vector<std::vector<uint32_t>>> pattern2rightFrontiers_vec;

    std::vector<std::unordered_set<uint32_t>> highDegreeNode;
    std::vector<std::vector<uint32_t>> highDegree_vec;

    std::unordered_map<uint32_t, int> lastEdge2SzLim;

    std::vector<uint64_t> effectiveNeighborsForFunFact;

    double smallest_one_hop = 0; //the largest possible score within 1-hop of the context nodes
    double scoreAmountExceptContextNodes = 0;

    ///for result summarization. This can also save half of the time
    std::unordered_map<size_t, int> resSz2resIdx;
    std::vector<std::vector<result_struct>> resIdx2resStruct;
    std::vector<std::vector<result_struct>> repeated_resIdx2resStruct;

public:
    /** Basic operations of priority queue. **/
    bool insert_topk_relevant_pattern(result_struct &r) {
        if (topk_relevant_pq.size() < top_pattern_k) {
            topk_relevant_pq.push(r);
            return true;
        } else if (topk_relevant_pq.top().patternRelevance < r.patternRelevance) {
            topk_relevant_pq.pop();
            topk_relevant_pq.push(r);
            return true;
        }
        return false;
    }

    /** This will be used for pruning. We need to look-up for the worst score. **/
    double current_kth_relevance_score() {
        if (topk_relevant_pq.size() < top_pattern_k) {
            return -1;
        } else {
            return topk_relevant_pq.top().patternRelevance;
        }
    }

public:

    /**Enumerate the subspace of the entity of interest.*/
    double
    enumSubspace(GraphManager &gm,
                 uint32_t target,
                 std::unordered_set<uint32_t> &peers,
                 uint32_t current_node, double patternRelevance,
                 std::string &patternSig);

    /**Calculate the exceptional facts around an entity
     * The input the attribute_value pairs.
     * It can be combinations of more than on attribute_value pair*/
    uint32_t mostFrequentValue = 0;
    double freqMFV = 0;
    int currentCount = 0;
    double
    exceptionalScoreCal(GraphManager &gm,
                        uint32_t target,
                        std::unordered_set<uint32_t> &peerEntities,
                        uint32_t attributeEid,
                        uint32_t attributeValueNode);


    void khopReverseExpandingFromContext(GraphManager &gm, uint32_t target, int hop_num);

/** Baseline solution.*/
    void enumPathPatternWithReverseExpansion_naive(GraphManager &gm, uint32_t target, int hop_num = 2);


    void backWardFromNode_Naive(GraphManager &gm, std::vector<uint32_t> &tmpEdgeTypes,
                                std::unordered_set<std::string> &visitedEdgeSeq,
                                int node2add, int target, double path_score);

    void emitResult(GraphManager &gm, std::string &edgeSeqSig,
                    std::unordered_set<uint32_t> &peers, uint32_t target,
                    uint32_t contextNid, double pScore);

    void emitResultJaccardSim(GraphManager &gm, std::string &edgeSeqSig,
                              std::unordered_set<uint32_t> &_fset,
                              uint32_t target,
                              uint32_t contextNid,
                              double pScore, double JDiffLim = 0.1);

/** Use backward optimization.*/
public:
    bool naiveDFS = false;
    std::vector<std::string> nodeSeqLabels_DFS;
    //DFS method
    std::unordered_map<uint32_t, std::vector<uint64_t>> hdNid2EffectNeighbors;

    void backWardFromNode_DFS_Reverse(GraphManager &gm,
                                      std::vector<uint32_t> &tmpEdgeTypes,
                                      std::unordered_set<std::string> &visitedEdgeSeq,
                                      double sz_exp,
                                      int node2add, int target, int hop_num,
                                      double path_score);

    double pathPatternFromContextNode_DFS(GraphManager &gm,
                                          std::vector<uint32_t> &reversedEdgeSeq,
                                          std::vector<std::unordered_set<uint32_t>> &levelFrontiers,
                                          int mid_idx, int target, int source);

    bool DFS_pattern_score(GraphManager &gm,
                           std::vector<uint32_t> &reversedEdgeSeq,
                           int edge_ptr,
                           std::vector<std::unordered_set<uint32_t>> &levelFrontiers,
                           std::vector<double> &levelScore,
                           std::vector<double> &numFrontier,
                           uint32_t cf, uint32_t source, bool &discard);

    bool testFrontierSzLimit(GraphManager &gm, int peerSz, uint32_t last_edge);

    //add pattern reachability
    void khopReverseExpandingFromContext_StrOnly(GraphManager &gm, uint32_t target, int hop_num);

    void enumPathPatternWithReverseExpansion_DFS_Reverse(GraphManager &gm, uint32_t target, int hop_num = 2);


    /** Competitors*/
public:
    void oneHopNeighborFind(GraphManager &gm, uint32_t target);

    void emitResult_Blind(GraphManager &gm, std::string &edgeSeqSig,
                          std::unordered_set<uint32_t> &_fset,
                          uint32_t target,
                          uint32_t contextNid,
                          double pScore);

    double enumSubspace_blind(GraphManager
                              &gm,
                              uint32_t target,
                              std::unordered_set<uint32_t>
                              &peers,
                              uint32_t current_node,
                              std::string
                              &patternSig);

    void enumPathPatternBlind(GraphManager &gm, uint32_t target);

    void backWardFromNode_Blind_EndingNode(GraphManager &gm, std::vector<uint32_t> &tmpEdgeTypes,
                                std::unordered_set<std::string> &visitedEdgeSeq,
                                int node2add, int target);

    void backWardFromNode_Blind_EndingType(GraphManager &gm, std::vector<uint32_t> &tmpEdgeTypes,
                                           std::unordered_set<std::string> &visitedEdgeSeq,
                                           int node2add, int target);


/** Remember high degree nodes for a certain prefix.*/
//public:
//    void enumPathPatternWithReverseExpansion_DFS_Reverse_MemorizeHD(GraphManager &gm, uint32_t target, int hop_num = 2);
//
//    bool DFS_pattern_score_MemorizeHD(GraphManager &gm,
//                                      std::vector<uint32_t> &reversedEdgeSeq,
//                                      int edge_ptr,
//                                      std::vector<std::unordered_set<uint32_t>> &levelFrontiers,
//                                      std::vector<std::unordered_set<uint32_t>> &levelHDNodes,
//                                      std::vector<double> &levelScore,
//                                      std::vector<double> &numFrontier,
//                                      uint32_t cf, uint32_t source,
//                                      bool &discard,
//                                      bool parentHD, uint32_t parentHDId);
//
//    void backWardFromNode_DFS_Reverse_MemorizeHD(GraphManager &gm,
//                                                 std::vector<uint32_t> &tmpEdgeTypes,
//                                                 std::unordered_set<std::string> &visitedEdgeSeq,
//                                                 double sz_exp,
//                                                 int node2add, int target, int hop_num);
//
//    double pathPatternFromContextNode_DFS_MemorizeHD(GraphManager &gm,
//                                                     std::vector<uint32_t> &reversedEdgeSeq,
//                                                     std::vector<std::unordered_set<uint32_t>> &levelFrontiers,
//                                                     int mid_idx, int target, int source);


    static inline
    std::string pathPattern2SigWNodes(std::vector<std::string> &nodeTypes, std::vector<uint32_t> &edgeIds) {
        std::string sig;
        for (int i = 0; i < edgeIds.size(); ++i) {
            sig += nodeTypes[i] + "," + std::to_string(edgeIds[i]) + ",";
        }
        return sig;
    }

    static inline
    void printResQueue(GraphManager &gm, uint32_t target,
                       Res_PQ resPq) {
        std::vector<result_struct> final_res;
        while (!resPq.empty()) {
            final_res.push_back(resPq.top());
            resPq.pop();
        }
        int count = 0;
        reverse(final_res.begin(), final_res.end());
        for (auto &r : final_res) {
            count++;
            //print Path Pattern
            auto eSeq = gm.sig2pathPattern(r.pattern_sig);

            std::string printPattern = r.nodeSeqLabels.front();

            for (int i = 0; i < eSeq.size(); ++i) {
                // edge
                if (eSeq[i] % 2 == 0)
                    printPattern.append("-[" + gm.typeId2Name[eSeq[i]] + "]->");
                else
                    printPattern.append("<-[" + gm.typeId2Name[gm.reverseEdgeTypeId(eSeq[i])] + "]-");
                // node
                printPattern.append(r.nodeSeqLabels[i + 1]);
            }


            printf("Path Pattern = %s, Entity Size = %ld, Relevance Score = %.9lf,\n\tFact = <%s, %s, %s>, Current Count = %d, Exceptional Score = %.6lf, most_value = <%s>, freq = %lf.\n",
                   printPattern.c_str(),
                   r.peer_entities.size(),
                   r.patternRelevance,
                   gm.nid2IRI[target].c_str(), gm.typeId2Name[r.fact_edge].c_str(), gm.nid2IRI[r.fact_nodeVal].c_str(),
                   r.current_count,
                   r.excepScore,
                   gm.nid2IRI[r.mostFrequentValue].c_str(), r.freqMFV);
        }
    }


    inline uint32_t getReversePathType(GraphManager &gm, std::vector<uint32_t> &left2mid,
                                       const std::string &mid2right, std::vector<uint32_t> &reverseEdge) {
        std::istringstream iss(mid2right);
        std::vector<uint32_t> tmp;
        uint32_t token;
        while (iss >> token) {
            tmp.push_back(token);
        }
        uint32_t contextNid;
        for (int i = 0; i < tmp.size(); ++i) {
            auto eid = tmp[tmp.size() - i - 1];
            if (i == 0) {
                contextNid = eid;
                continue;
            }
            reverseEdge.push_back(gm.reverseEdgeTypeId(eid));
        }
        for (int i = 0; i < left2mid.size(); ++i) {
            auto eid = left2mid[left2mid.size() - i - 1];
            reverseEdge.push_back(gm.reverseEdgeTypeId(eid));
        }
        return contextNid;
    }

    inline std::string reverseEdgeSeqSig(GraphManager &gm, std::vector<uint32_t> &edgeSeq) {
        std::string res;
        for (int i = 0; i < edgeSeq.size(); ++i) {
            auto eid = edgeSeq[edgeSeq.size() - i - 1];
            auto reversed_eid = gm.reverseEdgeTypeId(eid);
            if (i != 0)
                res.append(" ");
            res.append(std::to_string(reversed_eid));
        }
        return res;
    }

    inline double jaccardDistance(std::unordered_set<uint32_t> &s1, std::unordered_set<uint32_t> &s2) {
        double commonNum = 0;
        if (s1.size() > s2.size()) {
            for (auto s : s2) {
                if (s1.count(s))
                    commonNum++;
            }

        } else {
            for (auto s : s1) {
                if (s2.count(s)) {
                    commonNum++;
                }
            }
        }

        return 1 - commonNum / ((double) s1.size() + (double) s2.size() - commonNum);

    }

};


#endif //PATHPATTERNMINING_EXCEPTIONALFACTMINING_H
