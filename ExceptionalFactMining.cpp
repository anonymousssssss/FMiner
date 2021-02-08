//
// Created by yangyueji on 9/13/19.
//

#include "ExceptionalFactMining.h"
#include <iostream>
#include <queue>

using namespace std;


double ExceptionalFactMining::exceptionalScoreCal(GraphManager &gm,
                                                  uint32_t target,
                                                  std::unordered_set<uint32_t> &peerEntities,
                                                  uint32_t attributeEid,
                                                  uint32_t attributeValueNode) {

    //string targetAttValSig = attrValue2SigString(attributeEid, attributeValueNode);
    freqMFV = 0;
    int counter = 0;
    unordered_map<uint32_t, int> nodeVal2count;
    nodeVal2count[attributeValueNode] = 1;
    double totalCount = 1;
    for (auto nid : peerEntities) {
        if (nid == target) continue;

        //auto edgePos = gm.firstEdgePos(nid, attributeEid);
        auto _key = combine2u32(nid, attributeEid);
        if (!gm.snidEid2pos.count(_key)) continue;
        auto edgePos = gm.snidEid2pos[_key];
        if (edgePos == -1) {
            //no such value, we ignore it. Because we only care about the frequency of the facts, not the peer entities.
            continue;
        }
//        auto expected_hit = combine2u32(attributeEid, attributeValueNode);
//        bool has_end_node = std::binary_search(gm.edges + edgePos, gm.edges + gm.nodes[nid+1], expected_hit);
//        uint32_t node_val;
//        if (has_end_node) {
//            node_val = attributeValueNode;
//        }
        bool has_end_node = false;
        uint32_t node_val;

        //only consider the end_node with the highest degree
//        auto end_node = extractLow32bits(gm.edges[edgePos]);
//        if (nodeVal2count.count(end_node)) {
//            nodeVal2count[end_node]++;
//        } else {
//            nodeVal2count[end_node] = 1;
//        }
//        totalCount++;
        unordered_set<uint32_t> node_vals;
        for (int i = edgePos; i < gm.nodes[nid + 1] && !has_end_node; ++i) {
            auto end_node = extractLow32bits(gm.edges[i]);
            auto eid = extractHigh32bits(gm.edges[i]);
            if (eid != attributeEid) break;
            if (end_node == attributeValueNode) {
                has_end_node = true;
                node_val = end_node;
                break;
            }

            if (i == edgePos)
                node_val = end_node;

//            if (i != edgePos) break;
        }


        if (nodeVal2count.count(node_val)) {
            nodeVal2count[node_val]++;
        } else {
            nodeVal2count[node_val] = 1;
        }
        totalCount++;

    }
    if (totalCount < this->num_peerEntity) return 0;
    int targetCount = nodeVal2count[attributeValueNode];
    double higherCountTotal = 0;
//    double tmp_higherCountTotal = 0;
    for (auto &kv_pair : nodeVal2count) {
        if (kv_pair.second > targetCount) {
            higherCountTotal += kv_pair.second;
//            tmp_higherCountTotal += kv_pair.second;
//            higherCountTotal += kv_pair.second * (kv_pair.second - targetCount);
        }

        if (kv_pair.second > freqMFV) {
            freqMFV = kv_pair.second;
            mostFrequentValue = kv_pair.first;
        }

    }
    currentCount = nodeVal2count[attributeValueNode];
    double exceptionScore = higherCountTotal / totalCount;
//    double exceptionScore = higherCountTotal / (totalCount * totalCount);


//    if (peerEntities.size() == 3039763) {
//        printf("Here = %lf, %lf, %lf\n", tmp_higherCountTotal, totalCount, tmp_higherCountTotal / totalCount);
//        printf("test again: %lf, exceptScore = %lf\n", higherCountTotal, exceptionScore);
//    }


    return exceptionScore;
//    return tmp_exceptionScore;
}

double ExceptionalFactMining::enumSubspace_blind(GraphManager &gm,
                                                 uint32_t target,
                                                 std::unordered_set<uint32_t> &peers,
                                                 uint32_t current_node,
                                                 std::string &patternSig) {
    int edge_start = gm.nodes[target];
    int adj_sz = gm.nodes[target + 1] - edge_start;
    double largestScore = 0;
    uint32_t largestEid = 0, largestNodeVal = 0;
    uint32_t largest_mostFrequentValue = 0;
    double largest_freqMFV = 0;

    for (int i = 0; i < adj_sz; ++i) {
        uint64_t edge = gm.edges[edge_start + i];
        auto eid = extractHigh32bits(edge);
        auto end_node = extractLow32bits(edge);

        ///Check if the edge is effective
        double possibleNidSz = gm.typeId2Count[eid].size();
        double ratio = possibleNidSz / peers.size();
        if (ratio < 0.9) {
            continue;
        }

        double exceptionScore = exceptionalScoreCal(gm, target, peers, eid, end_node);

        if (exceptionScore == 0) continue;
        //cout << exceptionScore << endl;

        if (res_pq.size() < global_topk) {
            result_struct _rs(exceptionScore, eid, end_node, current_node, patternSig, 1, peers);
            _rs.set_node_labels(gm, target, nodeSeqLabels_DFS);
            _rs.mostFrequentValue = mostFrequentValue;
            _rs.freqMFV = freqMFV;
            _rs.current_count = currentCount;
            res_pq.push(_rs);
        } else {
            auto &top_r = res_pq.top();
            if (top_r.excepScore < exceptionScore) {
                res_pq.pop();
                result_struct _rs(exceptionScore, eid, end_node, current_node, patternSig, 1, peers);
                _rs.set_node_labels(gm, target, nodeSeqLabels_DFS);
                _rs.mostFrequentValue = mostFrequentValue;
                _rs.freqMFV = freqMFV;
                _rs.current_count = currentCount;
                res_pq.push(_rs);
            }
        }
        if (exceptionScore > largestScore) {
            largestScore = exceptionScore;
            largestEid = eid;
            largestNodeVal = end_node;
            largest_mostFrequentValue = mostFrequentValue;
            largest_freqMFV = freqMFV;
        }
//        omp_unset_lock(&lock);
    }

    if (largestScore == 0) return 0;
    if (diverseTopk.size() < global_topk) {
        result_struct _rs(largestScore, largestEid, largestNodeVal, current_node, patternSig, 1, peers);
        _rs.set_node_labels(gm, target, nodeSeqLabels_DFS);
        _rs.mostFrequentValue = largest_mostFrequentValue;
        _rs.freqMFV = largest_freqMFV;
        diverseTopk.push(_rs);
    } else {
        auto &top_r = diverseTopk.top();
        if (top_r.excepScore < largestScore) {
            diverseTopk.pop();
            result_struct _rs(largestScore, largestEid, largestNodeVal, current_node, patternSig, 1,
                              peers);
            _rs.set_node_labels(gm, target, nodeSeqLabels_DFS);
            _rs.mostFrequentValue = largest_mostFrequentValue;
            _rs.freqMFV = largest_freqMFV;
            diverseTopk.push(_rs);

        }
    }
    return 1;
}


double ExceptionalFactMining::enumSubspace(GraphManager &gm,
                                           uint32_t target,
                                           std::unordered_set<uint32_t> &peers,
                                           uint32_t current_node, double patternRelevance,
                                           std::string &patternSig) {


//    priority_queue<scoreFact, vector<scoreFact>, greater<scoreFact>> topk_pq;

    int edge_start = gm.nodes[target];
    int adj_sz = gm.nodes[target + 1] - edge_start;
    double largestScore = 0;
    uint32_t largestEid = 0, largestNodeVal = 0;
    uint32_t largest_mostFrequentValue = 0;
    double largest_freqMFV = 0;

//    omp_lock_t lock;
//    omp_init_lock(&lock);
    bool handled = !effectiveNeighborsForFunFact.empty();
//bool handled = false;
//#pragma omp for schedule(static)
    for (int i = 0; i < adj_sz; ++i) {
        uint64_t edge;
        if (!handled)
            edge = gm.edges[edge_start + i];
        else {
            if (i >= effectiveNeighborsForFunFact.size()) break;
            edge = effectiveNeighborsForFunFact[i];
        }

        auto eid = extractHigh32bits(edge);
        auto end_node = extractLow32bits(edge);

        ///check if the end_node is significant globally
        // commented for top-k search
//        if (!nodeScore.count(end_node) || nodeScore[end_node] < relevancy_limit) continue;
        if (!nodeScore.count(end_node) || nodeScore[end_node] < current_kth_relevance_score()) continue;


        if (!handled)
            effectiveNeighborsForFunFact.push_back(edge);

        ///Check if the edge is effective
        double possibleNidSz = gm.typeId2Count[eid].size();
        double ratio = possibleNidSz / peers.size();
        if (ratio < 0.8) {
            //cout << "Poss = "<< possibleNidSz << ", ratio = " << ratio << ", !!!!!!!!!!!!!!!" << endl;
            continue;
        }

        double exceptionScore = exceptionalScoreCal(gm, target, peers, eid, end_node);

        if (exceptionScore == 0) continue;
        //cout << exceptionScore << endl;
//        string res = gm.typeId2Name[eid] + ": " + gm.nid2IRI[end_node];
//        if (topk_pq.size() < topk_exceptionalFact) {
//            topk_pq.push(make_pair(exceptionScore, res));
//        } else {
//            auto &topPair = topk_pq.top();
//            if (topPair.first < exceptionScore) {
//                topk_pq.pop();
//                topk_pq.push(make_pair(exceptionScore, res));
//            }
//        }
//        result_struct _rs(exceptionScore, eid, end_node, current_node, patternSig, patternRelevance, peers);

//        omp_set_lock(&lock);
        if (res_pq.size() < global_topk) {
            result_struct _rs(exceptionScore, eid, end_node, current_node, patternSig, patternRelevance, peers);
            _rs.set_node_labels(gm, target, nodeSeqLabels_DFS);
            _rs.mostFrequentValue = mostFrequentValue;
            _rs.freqMFV = freqMFV;
            res_pq.push(_rs);
        } else {
            auto &top_r = res_pq.top();
            if (top_r.excepScore < exceptionScore) {
                res_pq.pop();
                result_struct _rs(exceptionScore, eid, end_node, current_node, patternSig, patternRelevance, peers);
                _rs.set_node_labels(gm, target, nodeSeqLabels_DFS);
                _rs.mostFrequentValue = mostFrequentValue;
                _rs.freqMFV = freqMFV;

                res_pq.push(_rs);
            }
        }
        if (exceptionScore > largestScore) {
            largestScore = exceptionScore;
            largestEid = eid;
            largestNodeVal = end_node;
            largest_mostFrequentValue = mostFrequentValue;
            largest_freqMFV = freqMFV;
        }
//        omp_unset_lock(&lock);
    }

    if (largestScore == 0) return 0;
    if (diverseTopk.size() < global_topk) {
        result_struct _rs(largestScore, largestEid, largestNodeVal, current_node, patternSig, patternRelevance, peers);
        _rs.set_node_labels(gm, target, nodeSeqLabels_DFS);
        _rs.mostFrequentValue = largest_mostFrequentValue;
        _rs.freqMFV = largest_freqMFV;
        diverseTopk.push(_rs);
    } else {
        auto &top_r = diverseTopk.top();
        if (top_r.excepScore < largestScore) {
            diverseTopk.pop();
            result_struct _rs(largestScore, largestEid, largestNodeVal, current_node, patternSig, patternRelevance,
                              peers);
            _rs.set_node_labels(gm, target, nodeSeqLabels_DFS);
            _rs.mostFrequentValue = largest_mostFrequentValue;
            _rs.freqMFV = largest_freqMFV;
            diverseTopk.push(_rs);

        }
    }
    return 1;

}

void ExceptionalFactMining::emitResult_Blind(GraphManager &gm, std::string &edgeSeqSig,
                                             std::unordered_set<uint32_t> &_fset,
                                             uint32_t target,
                                             uint32_t contextNid,
                                             double pScore) {
    //auto tmp_sig = gm.sig2labelSeq(edgeSeqSig);

    result_struct _tmp_rs(0, 0, 0, contextNid, edgeSeqSig, pScore);
//                    if (resSz2resIdx.count(_fset.size())) {
    int res_idx;
    if (resSz2resIdx.count(_fset.size())) {
        res_idx = resSz2resIdx[_fset.size()];
    } else {
        res_idx = resIdx2resStruct.size();
        resSz2resIdx[_fset.size()] = res_idx;
        resIdx2resStruct.resize(res_idx + 1);
        repeated_resIdx2resStruct.resize(res_idx + 1);
    }

    bool summarized = false;
    for (int i = 0; i < resIdx2resStruct[res_idx].size(); ++i) {
        if (_fset == resIdx2resStruct[res_idx][i].peer_entities) {
            _tmp_rs.fact_nodeVal = i;
            summarized = true;
        }
    }

    if (summarized) {
        //resIdx2resStruct[res_idx].push_back(_tmp_rs);
        repeated_resIdx2resStruct[res_idx].push_back(_tmp_rs);
        return;
    }

    //no repeat
    _tmp_rs.fact_nodeVal = resIdx2resStruct[res_idx].size();
    _tmp_rs.peer_entities = _fset;
    resIdx2resStruct[res_idx].push_back(_tmp_rs);
    repeated_resIdx2resStruct[res_idx].push_back(_tmp_rs);
//                    }

    enumSubspace_blind(gm, target, _fset, contextNid,edgeSeqSig);
}


void ExceptionalFactMining::emitResult(GraphManager &gm, std::string &edgeSeqSig,
                                       std::unordered_set<uint32_t> &_fset,
                                       uint32_t target,
                                       uint32_t contextNid,
                                       double pScore) {

    result_struct _tmp_rs(0, 0, 0, contextNid, edgeSeqSig, pScore);
//                    if (resSz2resIdx.count(_fset.size())) {

    // insert to top-k pattern relevance queue before fed into OF finding
    if ( !insert_topk_relevant_pattern(_tmp_rs) )
        return;

//    insert_topk_relevant_pattern(_tmp_rs);
//    return;

    int res_idx;
    if (resSz2resIdx.count(_fset.size())) {
        res_idx = resSz2resIdx[_fset.size()];
    } else {
        res_idx = resIdx2resStruct.size();
        resSz2resIdx[_fset.size()] = res_idx;
        resIdx2resStruct.resize(res_idx + 1);
        repeated_resIdx2resStruct.resize(res_idx + 1);
    }

    bool summarized = false;
    for (int i = 0; i < resIdx2resStruct[res_idx].size(); ++i) {
        if (_fset == resIdx2resStruct[res_idx][i].peer_entities) {
            _tmp_rs.fact_nodeVal = i;
            summarized = true;
        }
    }

    if (summarized) {
        //resIdx2resStruct[res_idx].push_back(_tmp_rs);
        repeated_resIdx2resStruct[res_idx].push_back(_tmp_rs);
        return;
    }

    //no repeat
    _tmp_rs.fact_nodeVal = resIdx2resStruct[res_idx].size();
    _tmp_rs.peer_entities = _fset;
    resIdx2resStruct[res_idx].push_back(_tmp_rs);
    repeated_resIdx2resStruct[res_idx].push_back(_tmp_rs);

    enumSubspace(gm, target, _fset, contextNid,
                 pScore, edgeSeqSig);
}


void ExceptionalFactMining::enumPathPatternWithReverseExpansion_naive(GraphManager &gm, uint32_t target, int hop_num) {
    // uncomment below to disable VR optimization
//    naiveDFS = true; // moved out

    uint64_t inner_s = getTime();
    double largestScore = 0;
    for (auto kv_pair : nodeScore) {
        if (!contextNodes.count(kv_pair.first)) {
            scoreAmountExceptContextNodes += kv_pair.second;
        } else {
            if (kv_pair.second > largestScore)
                largestScore = kv_pair.second;
        }

    }

    //khopReverseExpandingFromContext(gm, target, hop_num);

    unordered_map<size_t, unordered_map<uint32_t, string>> peerSz2firstPthIns;

    /// The string refers to the path instance.
    /// From path instance, we can infer the corresponding path pattern.
    /// The second element is the signature of path instance.
    typedef pair<double, pair<double, string>> patternSeqScore;
    /// By default, this is a max heap that is needed.
    priority_queue<patternSeqScore, vector<patternSeqScore>> pathQueue;
    vector<uint32_t> initPath;
    initPath.push_back(target);
    string _sig = gm.pathPattern2sig(initPath);
    double _initScore = 1;
    pathQueue.push(make_pair(_initScore, make_pair(1, _sig)));

    /// For repetition checking
    unordered_map<uint32_t, unordered_set<string>> endNode2pathSig;
    unordered_set<string> visitedEdgeSeq;

    /// For preventing the path patterns leading to unchanged peer entities
    unordered_map<uint32_t, unordered_map<string, unordered_set<uint32_t>>> node2pattern2peerEntities;


    /// need to handle the single target node without expansion
//    vector<uint32_t> emptyEdgeTypes;
//    auto test1 = getTime();
//    backWardFromNode_Naive(gm, emptyEdgeTypes, visitedEdgeSeq, 0, target, target, hop_num);
//    auto test2 = getTime();
//    backwardFromNodeTime += getInterval(test1, test2);

    while (!pathQueue.empty()) {
        auto tLimit = getTime();
        if (getInterval(inner_s, tLimit) >= TimeLimit) {
            InnerProcessingTime = TimeLimit;
            return;
        }

        auto curr_pair = pathQueue.top();
        pathQueue.pop();

        /**The path score decreases monotonicity.**/
        if (curr_pair.first < current_kth_relevance_score()) {
            // This means the operation
            break;
        }


        auto pathInsSig = curr_pair.second;


        ///Get path patterns
        auto pathIns = gm.sig2pathPattern(pathInsSig.second);
        vector<string> nodeTypePattern;
        vector<uint32_t> edgeTypePattern;
        uint32_t current_node = pathIns.back();

        for (int i = 0; i < pathIns.size(); ++i) {
            if (i % 2 == 0) {
                //it is a node
                auto nodeTypes = gm.nid2types[pathIns[i]]; // node type
                // if this is the case, then the node has no type. We fix this node
                if (nodeTypes.empty())
                    nodeTypePattern.emplace_back("ANY");
                for (auto &np : nodeTypes) {
                    nodeTypePattern.push_back(np);
                    break;
                }
            } else {
                edgeTypePattern.push_back(pathIns[i]);
            }
        }

        ///Expand the node
        size_t edge_start = gm.nodes[current_node];
        size_t adj_sz = gm.nodes[current_node + 1] - edge_start;

        for (int i = 0; i < adj_sz; ++i) {
            auto edge = gm.edges[edge_start + i];
            uint32_t eid = extractHigh32bits(edge);
            uint32_t end_node = extractLow32bits(edge);
            double sz_exp = nodeTypePattern.size();

            if (sz_exp >= path_length_limit) continue;
            if (end_node == target) continue; //inside node cannot be either contextNode or target node
            if (current_node == end_node) continue; //no self-cycle
            if (!nodeScore.count(end_node)) {
                continue; // the node score is too low. cancel path instance level relevancy limit
            }
            double end_node_score = nodeScore[end_node] * pow(base_e, -sz_exp);
            double path_score = end_node_score < curr_pair.first ? end_node_score : curr_pair.first;

            if (path_score <= current_kth_relevance_score())
                continue;

            // commented for top-k search
//            if (path_score < relevancy_limit) {
//                continue; // cancel path instance level relevancy
//            }

            ///Check possible repetitions
            vector<uint32_t> tmpEdgeTypes = edgeTypePattern;
            tmpEdgeTypes.push_back(eid);

            auto sigRep = pathPattern2SigWNodes(nodeTypePattern, tmpEdgeTypes);


            if (endNode2pathSig.count(end_node) && endNode2pathSig[end_node].count(sigRep))
                continue;
            endNode2pathSig[end_node].insert(sigRep);


            auto test1 = getTime();
            backWardFromNode_Naive(gm, tmpEdgeTypes, visitedEdgeSeq,
                    end_node, target, path_score);
            auto test2 = getTime();
            backwardFromNodeTime += getInterval(test1, test2);

            if (contextNodes.count(end_node)) continue;

            double possible_largest = largestScore * pow(base_e, -(sz_exp + 1));
            possible_largest = path_score < possible_largest ? path_score : possible_largest;
//            if (possible_largest < relevancy_limit) {
//  //              continue; //cancel path instance level relevancy limit
//            }

            vector<uint32_t> tmpIns = pathIns;
            tmpIns.push_back(eid);
            tmpIns.push_back(end_node);
            pathQueue.push(make_pair(path_score, make_pair(possible_largest, gm.pathPattern2sig(tmpIns))));
        }
    }


    uint64_t inner_e = getTime();

    ///print final info
//    cout << "\n#### Output global top-k exceptional facts" << endl;
//    printResQueue(gm, target, res_pq);

    //cout << "\n#### Output global Diverse top-k exceptional facts" << endl;
    //printResQueue(gm, target, diverseTopk);

    //printf("Inner time = %.3lf ms.\n", getInterval(inner_s, inner_e));
    InnerProcessingTime = getInterval(inner_s, inner_e);

}


void ExceptionalFactMining::emitResultJaccardSim(GraphManager &gm, std::string &edgeSeqSig,
                                                 std::unordered_set<uint32_t> &_fset,
                                                 uint32_t target,
                                                 uint32_t contextNid,
                                                 double pScore, double JsimDis) {
    //auto tmp_sig = gm.sig2labelSeq(edgeSeqSig);

    result_struct _tmp_rs(0, 0, 0, contextNid, edgeSeqSig, pScore);
//                    if (resSz2resIdx.count(_fset.size())) {
    int res_idx;
    int lowSz = _fset.size() * (1 - JsimDis);
    lowSz = lowSz > 0 ? lowSz : 0;
    int highSz = _fset.size() * (1 + JsimDis);
    bool summarized = false;
    for (size_t i = lowSz; i <= highSz; ++i) {
        if (resSz2resIdx.count(i)) {
            res_idx = resSz2resIdx[i];

//        else {
//            res_idx = resIdx2resStruct.size();
//            resSz2resIdx[_fset.size()] = res_idx;
//            resIdx2resStruct.resize(res_idx + 1);
//            repeated_resIdx2resStruct.resize(res_idx + 1);
//        }


            for (int j = 0; j < resIdx2resStruct[res_idx].size(); ++j) {
                double diff = jaccardDistance(_fset, resIdx2resStruct[res_idx][j].peer_entities);
                if (diff < JsimDis) {
                    summarized = true;
                    break;
                }
            }

            if (summarized) {
                //resIdx2resStruct[res_idx].push_back(_tmp_rs);
                repeated_resIdx2resStruct[res_idx].push_back(_tmp_rs);
                return;
            }
        }
    }

    //not able to be summarized
    if (!resSz2resIdx.count(_fset.size())) {
        res_idx = resIdx2resStruct.size();
        resSz2resIdx[_fset.size()] = res_idx;
        resIdx2resStruct.resize(res_idx + 1);
        repeated_resIdx2resStruct.resize(res_idx + 1);
    } else {
        res_idx = resSz2resIdx[_fset.size()];
    }

    //no repeat
    _tmp_rs.fact_nodeVal = resIdx2resStruct[res_idx].size();
    _tmp_rs.peer_entities = _fset;
    resIdx2resStruct[res_idx].push_back(_tmp_rs);
    repeated_resIdx2resStruct[res_idx].push_back(_tmp_rs);
//                    }

    enumSubspace(gm, target, _fset, contextNid,
                 pScore, edgeSeqSig);
}


void ExceptionalFactMining::khopReverseExpandingFromContext(GraphManager &gm, uint32_t target, int hop_num) {
    auto s_time = getTime();
    hopNid2Parent.clear();
    hopNid2Parent.resize(hop_num);
    highDegreeNode.resize(hop_num); // used for reachability
    unordered_set<uint32_t> frontiers = contextNodes;

    for (int i = 0; i < hop_num && !frontiers.empty(); ++i) {
        for (auto f : frontiers) {
            if (f == target) continue;//end node can be target
            int edge_start = gm.nodes[f];
            int adj_sz = gm.nodes[f + 1] - edge_start;
            if (adj_sz >= degree_threshold) {
                if (i == 0) {
                    //identify high degree context nodes
                    double _hd_score = nodeScore[f];
                    smallest_one_hop = smallest_one_hop > _hd_score ? smallest_one_hop : _hd_score;
                    highDegreeContextNodes.insert(f);
                }
                continue; // skip high degree node
            }


            //if (!nodeScore.count(f) || nodeScore[f] < relevancy_limit) continue;
            double f_score = 0;
            if (nodeScore.count(f)) f_score = nodeScore[f];

            for (int j = 0; j < adj_sz; ++j) {
                auto edge = gm.edges[edge_start + j];
                auto e_nid = extractLow32bits(edge);
                if (contextNodes.count(e_nid)) continue;
                if (f == e_nid) continue;
                if (gm.nodes[e_nid + 1] - gm.nodes[e_nid] >= degree_threshold)
                    highDegreeNode[i].insert(e_nid);
                //double c_score = nodeScore[e_nid];
                auto origin_eType = extractHigh32bits(edge);
                auto e_type = gm.reverseEdgeTypeId(origin_eType);
//                auto reverse_edge = combine2u32(e_type, f);
                string etoken = to_string(e_type);
                //string reverse_etoken = to_string(extractHigh32bits(edge));
//                if (!nodeScore.count(e_nid) || nodeScore[e_nid] < relevancy_limit) continue;


                if (i == 0) {
                    string _etoken = etoken + " " + to_string(f);
                    double p_score = f_score / base_e;
                    //if (nodeScore.count(e_nid) && nodeScore[e_nid] >= relevancy_limit && p_score >= relevancy_limit) {
                    ///legal path to keep
                    double nid_score = 0;
                    if (nodeScore.count(e_nid)) nid_score = nodeScore[e_nid];
                    p_score = p_score < nid_score ? p_score : nid_score;
                    if (!hopNid2Parent[i][e_nid].count(_etoken) || hopNid2Parent[i][e_nid][_etoken] < p_score) {
                        hopNid2Parent[i][e_nid][_etoken] = p_score;
                    }
//                    pattern2leftnode[_etoken].insert(e_nid);
                    //in the first hop. It is always correct.
                    smallest_one_hop = smallest_one_hop > p_score ? smallest_one_hop : p_score;
                    //}

                    if (!pattern2rightFrontiers_vec.count(_etoken)) {
                        unordered_set<uint32_t> tmp;
                        int _epos = gm.firstEdgePos(f, origin_eType);
                        if (_epos == -1) continue;
                        for (int k = _epos; k < gm.nodes[f + 1]; ++k) {
                            auto _tmp_edge = gm.edges[k];
                            if (extractHigh32bits(_tmp_edge) != origin_eType) break;
//                            pattern2rightFrontiers[_etoken].back().insert(extractLow32bits(_tmp_edge));
                            tmp.insert(extractLow32bits(_tmp_edge));
                        }
                        vector<uint32_t> t_vec;
                        pattern2rightFrontiers_vec[_etoken].push_back(t_vec);
                        for (auto t : tmp) {
                            pattern2rightFrontiers_vec[_etoken].back().push_back(t);
                        }
                    }

//                    if (pattern2rightFrontiers.count(_etoken)) {
//                        pattern2rightFrontiers[_etoken].back().insert(e_nid);
//                    } else {
//                        unordered_set<uint32_t> tmp;
//                        tmp.insert(e_nid);
////                        pattern2rightFrontiers[_etoken] = pattern2rightFrontiers[etoken];
//                        pattern2rightFrontiers[_etoken].push_back(tmp);
//                    }
                } else {
                    for (auto &s_pair : hopNid2Parent[i - 1][f]) {
                        string _etoken = etoken + " " + s_pair.first;
                        double p_score = s_pair.second / base_e;
//                        if (nodeScore.count(e_nid) && nodeScore[e_nid] >= relevancy_limit &&
//                            p_score >= relevancy_limit) {
                        ///legal path to keep
                        double nid_score = 0;
                        if (nodeScore.count(e_nid)) nid_score = nodeScore[e_nid];
                        p_score = p_score < nid_score ? p_score : nid_score;
                        if (!hopNid2Parent[i][e_nid].count(_etoken) || hopNid2Parent[i][e_nid][_etoken] < p_score) {
                            hopNid2Parent[i][e_nid][_etoken] = p_score;
                        }
                        //}
//                        pattern2leftnode[_etoken].insert(e_nid);
//                        if (!pattern2rightFrontiers.count(_etoken)) {
//                            //recover last layer frontiers.
//
//                        }
                        if (!pattern2rightFrontiers_vec.count(_etoken)) {
                            /// find the level frontiers including highdegree nodes.
                            pattern2rightFrontiers_vec[_etoken] = pattern2rightFrontiers_vec[s_pair.first];
//                            if (pattern2rightFrontiers_vec[s_pair.first].size() == 1 &&
//                                pattern2rightFrontiers_vec[s_pair.first].back()[0] == target) {
//                                continue; // do not include target
//                            }
//
//                            for (int l = 0; l < pattern2rightFrontiers_vec[s_pair.first].size(); ++l) {
//                                if (l != pattern2rightFrontiers_vec[s_pair.first].size() - 1)
//                                    pattern2rightFrontiers_vec[_etoken].push_back(
//                                            pattern2rightFrontiers_vec[s_pair.first][l]);
//                                else {
//                                    vector<uint32_t> t_vec;
//                                    pattern2rightFrontiers_vec[_etoken].push_back(t_vec);
//                                    for (auto k : pattern2rightFrontiers_vec[s_pair.first][l]) {
//                                        if (k != target && !contextNodes.count(k))
//                                            pattern2rightFrontiers_vec[_etoken].back().push_back(k);
//                                    }
//                                }
//
//                            }
//                            pattern2rightFrontiers_vec[_etoken].back().erase(target);
                            unordered_set<uint32_t> tmp;
//                            auto & _toExpand = pattern2rightFrontiers[_etoken].back();
                            int p_idx = (int) pattern2rightFrontiers_vec[_etoken].size() - 1;
                            //pattern2rightFrontiers_vec[_etoken].push_back(tmp);
//                            if (gm.typeId2Name[e_type] == "native language") {
//                                cout << "here sz = " << pattern2rightFrontiers_vec[_etoken][p_idx].size() << endl;
//                            }
                            for (auto m : pattern2rightFrontiers_vec[_etoken][p_idx]) {
//                                if (_etoken == "54 55 3198") {
//                                    cout << m << ", "<<gm.nid2IRI[m] << endl;
//                                }
//                                if (gm.typeId2Name[e_type] == "native language") {
//                                    cout << m << ", " << gm.nid2IRI[m] << endl;
//                                }

                                int _epos = gm.firstEdgePos(m, origin_eType);
                                if (_epos == -1) continue;
                                for (int k = _epos; k < gm.nodes[m + 1]; ++k) {
                                    auto _tmp_edge = gm.edges[k];
                                    if (extractHigh32bits(_tmp_edge) != origin_eType) break;
                                    auto _end_node = extractLow32bits(_tmp_edge);
//                                    if (!contextNodes.count(_end_node))
                                    tmp.insert(_end_node);
                                }
                            }
                            vector<uint32_t> t_vec;
                            pattern2rightFrontiers_vec[_etoken].push_back(t_vec);
                            for (auto t : tmp) {
                                pattern2rightFrontiers_vec[_etoken].back().push_back(t);

                            }

                        }
                    }
                }
            }
        }


        frontiers.clear();
        for (const auto &nextf : hopNid2Parent[i])
            frontiers.insert(nextf.first);
    }
    auto e_time1 = getTime();
    hopNid2Parent_vec.resize(hopNid2Parent.size());
    for (int i = 0; i < hopNid2Parent.size(); ++i) {
        for (auto &kv_pair : hopNid2Parent[i]) {
            vector<pair<double, string>> tmp;
            for (auto &_kv : kv_pair.second) {
                tmp.emplace_back(_kv.second, _kv.first);
            }
            std::sort(tmp.begin(), tmp.end(), greater<pair<double, string>>());
            hopNid2Parent_vec[i][kv_pair.first] = tmp;
        }

    }

    highDegree_vec.resize(hop_num);
    for (int n = 0; n < highDegreeNode.size(); ++n) {
        for (auto nid : highDegreeNode[n]) {
            highDegree_vec[n].push_back(nid);
        }
    }

    //already done
//    for (auto &kv_pair : pattern2rightFrontiers) {
//        for (auto v : kv_pair.second) {
//            vector<uint32_t> tmp_vec(v.begin(), v.end());
//            pattern2rightFrontiers_vec[kv_pair.first].push_back(tmp_vec);
//        }
//    }

    auto e_time = getTime();
//    printf("reverse expansion time = %lf ms. \n", getInterval(s_time, e_time1));
//    cout << endl;
}


void ExceptionalFactMining::backWardFromNode_Naive(GraphManager &gm, std::vector<uint32_t> &tmpEdgeTypes,
                                                   std::unordered_set<std::string> &visitedEdgeSeq,
                                                   int node2add, int target, double path_score) {
    if (node2add == target || !contextNodes.count(node2add)) return;

    vector<unordered_set<uint32_t>> levelFrontiers;
    auto reverseEdgePattern = gm.reversePatternSequence(tmpEdgeTypes);
    auto edgeSeqSig = reverseEdgeSeqSig(gm, reverseEdgePattern);
    if (visitedEdgeSeq.count(edgeSeqSig)) {
        return;
    }
    visitedEdgeSeq.insert(edgeSeqSig);
    levelFrontiers.resize(reverseEdgePattern.size() + 1);
    levelFrontiers.back().insert(node2add);

//    double pScore = pathPatternFromContextNode(gm,
//                                               reverseEdgePattern,
//                                               levelFrontiers,
//                                               (int) tmpEdgeTypes.size(),
//                                               target,
//                                               node2add);

    //by DFS
    nodeSeqLabels_DFS.clear();
    double pScore = pathPatternFromContextNode_DFS(gm, reverseEdgePattern,
                                                   levelFrontiers,
                                                   (int) reverseEdgePattern.size(),
                                                   target, node2add);

    //added for top-k search
    path_score = pScore < path_score? pScore : path_score;

    if (pScore > 0 && levelFrontiers.front().size() >= num_peerEntity) {
        auto t1 = getTime();
        emitResult(gm, edgeSeqSig, levelFrontiers.front(), target, node2add, path_score);
        auto t2 = getTime();
        emitResTime += getInterval(t1, t2);
    }

}

bool ExceptionalFactMining::DFS_pattern_score(GraphManager &gm, std::vector<uint32_t> &reversedEdgeSeq, int edge_ptr,
                                              std::vector<std::unordered_set<uint32_t>> &levelFrontiers,
                                              std::vector<double> &levelScore,
                                              std::vector<double> &numFrontier,
                                              uint32_t cf, uint32_t source, bool &discard) {
    if (discard) return false;
    if (edge_ptr == reversedEdgeSeq.size()) {
        levelFrontiers[0].insert(cf); // last edge
        if (!contextNodes.count(cf)) {
            numFrontier[0]++;
            if (nodeScore.count(cf))
                levelScore[0] += nodeScore[cf];
        }
        //test limit num
        if (!naiveDFS)
            if (testFrontierSzLimit(gm, (int) numFrontier[0], reversedEdgeSeq.back())) {
                discard = true;
                return false;
            }
        return true;
    }
    int levelFrontier_idx = (int) reversedEdgeSeq.size() - edge_ptr - 1;
    int edge_pos = gm.firstEdgePos(cf, reversedEdgeSeq[edge_ptr]);
    if (edge_pos == -1) return false;
    bool inserted = false;
    for (int i = edge_pos; i < gm.nodes[cf + 1]; ++i) {
        auto eid = extractHigh32bits(gm.edges[i]);
        if (eid != reversedEdgeSeq[edge_ptr]) break;
        uint32_t end_node = extractLow32bits(gm.edges[i]);
        if (end_node == cf) continue;
        //one edge look forward
        if (!naiveDFS)
            if (edge_ptr < reversedEdgeSeq.size() - 1) {
                uint32_t oneMoreEdge = reversedEdgeSeq[edge_ptr + 1];
                if (!gm.typeId2Count[oneMoreEdge].count(end_node)) {
                    continue;
                }
            }


        bool valid = levelFrontiers[levelFrontier_idx].count(end_node);

        if (!valid)
            valid = DFS_pattern_score(gm, reversedEdgeSeq, edge_ptr + 1,
                                      levelFrontiers, levelScore, numFrontier, end_node, source, discard);


        if (discard) {
            return false;
        }
        if (valid) {
            if (!inserted) {
                levelFrontiers[levelFrontier_idx + 1].insert(cf);
                inserted = true;
                if (!contextNodes.count(cf)) {
                    //this ensures the score calculation wont involve the context node
                    numFrontier[levelFrontier_idx + 1]++;
                    if (nodeScore.count(cf))
                        levelScore[levelFrontier_idx + 1] += nodeScore[cf];
                }

                //later for printing labels
                if (edge_ptr != 0 && edge_ptr + 1 + nodeSeqLabels_DFS.size() == reversedEdgeSeq.size()) {
                    string _label;
                    gm.getOneLabelRandom(cf, _label);
                    nodeSeqLabels_DFS.push_back(_label);
                }

                if (!naiveDFS)
                    if (edge_ptr > 0)
                        if (testFrontierSzLimit(gm, (int) numFrontier[levelFrontier_idx + 1],
                                                reversedEdgeSeq[edge_ptr - 1])) {
                            discard = true;
                            return false;
                        }

            }
        }
    }
    return inserted;
}


bool ExceptionalFactMining::testFrontierSzLimit(GraphManager &gm, int peerSz, uint32_t last_edge) {

    // used for top-k search
    if (topk_relevant_pq.size() < top_pattern_k)
        return false;// use the smallest top-k value for pruning.
    double kth_smallest = topk_relevant_pq.top().patternRelevance;

    int sz_lim = 0;
    last_edge = gm.reverseEdgeTypeId(last_edge);
    if (lastEdge2SzLim.count(last_edge))
        sz_lim = lastEdge2SzLim[last_edge];
    else {
        const auto &startNodes = gm.typeId2Count[last_edge];
        if (nodeScore.size() > startNodes.size()) {
            double s_amount = 0;
            for (auto sn : startNodes)
                if (nodeScore.count(sn) && !contextNodes.count(sn))
                    s_amount += nodeScore[sn];
// commented to use top-k search
//            s_amount = s_amount / relevancy_limit + 1;

            s_amount = s_amount / kth_smallest + 1;
            sz_lim = (int) s_amount;
            lastEdge2SzLim[last_edge] = (int) sz_lim;
        } else {
            double s_amount = 0;
            for (auto &kv_pair : nodeScore)
                if (startNodes.count(kv_pair.first) && !contextNodes.count(kv_pair.first))
                    s_amount += nodeScore[kv_pair.first];
              // commented to use top-k search
//            s_amount = s_amount / relevancy_limit + 1;
            s_amount = s_amount / kth_smallest + 1;

            sz_lim = (int) s_amount;
            lastEdge2SzLim[last_edge] = (int) sz_lim;
        }
    }

    return peerSz >= sz_lim;
}


double ExceptionalFactMining::pathPatternFromContextNode_DFS(GraphManager &gm, std::vector<uint32_t> &reversedEdgeSeq,
                                                             std::vector<std::unordered_set<uint32_t>> &levelFrontiers,
                                                             int mid_idx, int target, int source) {

    vector<double> levelScore(reversedEdgeSeq.size() + 1);
    vector<double> numFrontier(reversedEdgeSeq.size() + 1);
    levelScore.back() += nodeScore[source];
    numFrontier.back()++;

    bool discard = false;
    int rest_sz = (int) reversedEdgeSeq.size() - mid_idx;
    int levelF_idx = mid_idx - 1;
    for (auto cf : levelFrontiers[levelF_idx + 1]) {
        bool valid = DFS_pattern_score(gm, reversedEdgeSeq, rest_sz,
                                       levelFrontiers, levelScore, numFrontier, cf, source, discard);
        if (discard) {
            return -1;
        }
        if (valid) {
            levelFrontiers[levelF_idx + 1].insert(cf);
        }
    }


    if (levelFrontiers.front().size() < num_peerEntity)
        return -1;
    double s = 1;
    for (int i = 0; i < levelScore.size() - 1; ++i) {
        if (numFrontier[i] == 0) return -1;
        double tmp_s = levelScore[i] / numFrontier[i];
        tmp_s *= pow(base_e, -i);
        s = s < tmp_s ? s : tmp_s;

        // commented to use top-k search
//        if (s < relevancy_limit) return -1;
        if (s < current_kth_relevance_score()) return -1;

    }
    nodeSeqLabels_DFS.push_back(gm.nid2IRI[source]);
    return s;
}

void ExceptionalFactMining::khopReverseExpandingFromContext_StrOnly(GraphManager &gm, uint32_t target, int hop_num) {

    auto s_time = getTime();
    hopNid2Parent.clear();
    hopNid2Parent.resize(hop_num);
    //highDegreeNode.resize(hop_num); // used for reachability
    unordered_set<uint32_t> frontiers = contextNodes;

    for (int i = 0; i < hop_num && !frontiers.empty(); ++i) {
        for (auto f : frontiers) {
            if (f == target) continue;//end node can be target
            int edge_start = gm.nodes[f];
            int adj_sz = gm.nodes[f + 1] - edge_start;
            if (i != 0 && adj_sz >= degree_threshold) {
                continue; // skip high degree node
            }

            // commented for topk search
//            if (!nodeScore.count(f) || nodeScore[f] < relevancy_limit) continue;
            if (!nodeScore.count(f)) continue;

            double f_score = 0;
            if (nodeScore.count(f)) f_score = nodeScore[f];

            for (int j = 0; j < adj_sz; ++j) {
                auto edge = gm.edges[edge_start + j];
                auto e_nid = extractLow32bits(edge);
                if (contextNodes.count(e_nid)) continue;
                if (f == e_nid) continue;

                auto origin_eType = extractHigh32bits(edge);
                auto e_type = gm.reverseEdgeTypeId(origin_eType);
                string etoken = to_string(e_type);

                if (i == 0) {
                    string _etoken = etoken + " " + to_string(f);
                    double p_score = f_score / base_e;
                    // commented for topk search
//                    if (nodeScore.count(e_nid) && nodeScore[e_nid] >= relevancy_limit && p_score >= relevancy_limit) {
                    if (nodeScore.count(e_nid)) {

                            ///legal path to keep
                        double nid_score = 0;
                        nid_score = nodeScore[e_nid];
                        p_score = p_score < nid_score ? p_score : nid_score;
                        if (!hopNid2Parent[i][e_nid].count(_etoken) || hopNid2Parent[i][e_nid][_etoken] < p_score) {
                            hopNid2Parent[i][e_nid][_etoken] = p_score;
                        }

                        smallest_one_hop = smallest_one_hop > p_score ? smallest_one_hop : p_score;
                    }
                } else {
                    for (auto &s_pair : hopNid2Parent[i - 1][f]) {
                        string _etoken = etoken + " " + s_pair.first;
                        double p_score = s_pair.second / base_e;
                        // commented for top-k search.
//                        if (nodeScore.count(e_nid) && nodeScore[e_nid] >= relevancy_limit &&
//                            p_score >= relevancy_limit)
                        if (nodeScore.count(e_nid)) {
                            ///legal path to keep
                            double nid_score = 0;
                            if (nodeScore.count(e_nid)) nid_score = nodeScore[e_nid];
                            p_score = p_score < nid_score ? p_score : nid_score;
                            if (!hopNid2Parent[i][e_nid].count(_etoken) || hopNid2Parent[i][e_nid][_etoken] < p_score) {
                                hopNid2Parent[i][e_nid][_etoken] = p_score;
                            }
                        }
                    }
                }
            }
        }

        frontiers.clear();
        for (const auto &nextf : hopNid2Parent[i])
            frontiers.insert(nextf.first);
    }
    auto e_time1 = getTime();
    hopNid2Parent_vec.resize(hopNid2Parent.size());
    for (int i = 0; i < hopNid2Parent.size(); ++i) {
        for (auto &kv_pair : hopNid2Parent[i]) {
            vector<pair<double, string>> tmp;
            for (auto &_kv : kv_pair.second) {
                tmp.emplace_back(_kv.second, _kv.first);
            }
            std::sort(tmp.begin(), tmp.end(), greater<pair<double, string>>());
            hopNid2Parent_vec[i][kv_pair.first] = tmp;
        }

    }

    reverseExpansionTime += getInterval(s_time, e_time1);
    //printf("reverse expansion time = %lf ms. \n", getInterval(s_time, e_time1));

}

void ExceptionalFactMining::enumPathPatternWithReverseExpansion_DFS_Reverse(GraphManager &gm, uint32_t target,
                                                                            int hop_num) {
    //comment to test reverse search without optimizing pattern evaluation.
   // naiveDFS = true;

    uint64_t inner_s = getTime();
    khopReverseExpandingFromContext_StrOnly(gm, target, hop_num);

    unordered_map<size_t, unordered_map<uint32_t, string>> peerSz2firstPthIns;

    /// most relevant node type: used for avoiding repetitive computation
    unordered_map<uint32_t, string> most_rel_node_type;
    unordered_map<string, double> type_rel_score;

    /// The string refers to the path instance.
    /// From path instance, we can infer the corresponding path pattern.
    /// The second element is the signature of path instance.
    typedef pair<double, pair<double, string>> patternSeqScore;
    /// By default, this is a max heap that is needed.
    priority_queue<patternSeqScore, vector<patternSeqScore>> pathQueue;
    vector<uint32_t> initPath;
    initPath.push_back(target);
    string _sig = gm.pathPattern2sig(initPath);
    double _initScore = 1;
    pathQueue.push(make_pair(_initScore, make_pair(1, _sig)));

    /// For repetition checking
    unordered_map<uint32_t, unordered_set<string>> endNode2pathSig;
    unordered_set<string> visitedEdgeSeq;

    /// For preventing the path patterns leading to unchanged peer entities
    unordered_map<uint32_t, unordered_map<string, unordered_set<uint32_t>>> node2pattern2peerEntities;


    /// need to handle the single target node without expansion
    vector<uint32_t> emptyEdgeTypes;
    auto test1 = getTime();
    backWardFromNode_DFS_Reverse(gm, emptyEdgeTypes, visitedEdgeSeq, 0,
            target, target, hop_num, 1);
//    backWardFromNode_Naive(gm, emptyEdgeTypes, visitedEdgeSeq, 0, target, target, hop_num);
    auto test2 = getTime();
    backwardFromNodeTime += getInterval(test1, test2);
    while (!pathQueue.empty()) {
        //add a time limit
        auto tLimit = getTime();
        if (getInterval(inner_s, tLimit) >= TimeLimit) {
            InnerProcessingTime = TimeLimit;
            return;
        }
        auto curr_pair = pathQueue.top();
        pathQueue.pop();

        //early stop by topk search
        if (curr_pair.first < current_kth_relevance_score()) {
            break;
        }

        auto pathInsSig = curr_pair.second;

        ///Get path patterns
        auto pathIns = gm.sig2pathPattern(pathInsSig.second);
        vector<string> nodeTypePattern;
        vector<uint32_t> edgeTypePattern;
        uint32_t current_node = pathIns.back();

        bool abandon = false;
        for (int i = 0; i < pathIns.size(); ++i) {
            if (i % 2 == 0) {
                if (most_rel_node_type.count(pathIns[i])) {
                    // Already calculated.
                    if (most_rel_node_type[pathIns[i]].empty()) {
                        abandon = true;

                    } else {
                        nodeTypePattern.push_back(most_rel_node_type[pathIns[i]]);
                    }
                    break;
                }

                const auto& nodeTypes = gm.nid2types[pathIns[i]];
                if (nodeTypes.empty()) {
                    abandon = true;
                    break;
                }

                // impossible for score to be larger than 1
                double max_type_rel_score = 0;
                string most_rel_type;
                for (const auto& t : nodeTypes) {
                    double tmp_score = 0;
                    if (type_rel_score.count(t)) {
                        tmp_score = type_rel_score[t];
                    } else {
                        // consider this node type
                        double accum_count = 0;
                        for (uint32_t nid : gm.type2nid[t]) {
                            if (nid != target) {
                                accum_count += 1;
                                if (nodeScore.count(nid)) {
                                    tmp_score += nodeScore[nid];
                                }
                            }
                        }
                        if (accum_count > 0)
                            tmp_score /= accum_count;
                        type_rel_score[t] = tmp_score;
                    }

                    if (tmp_score > max_type_rel_score) {
                        max_type_rel_score = tmp_score;
                        most_rel_type = t;
                    }
                }

                if (max_type_rel_score == 0 || most_rel_type.empty()) {
                    abandon = true;
                    most_rel_node_type[pathIns[i]] = "";
                    break;
                }
                nodeTypePattern.push_back(most_rel_type);
                most_rel_node_type[pathIns[i]] = most_rel_type;
                /* Switch to use most specific node types as above. */
//                //it is a node
//                auto nodeTypes = gm.nid2types[pathIns[i]]; // node type
//                // if this is the case, then the node has no type. We fix this node
//                if (nodeTypes.empty())
//                    nodeTypePattern.emplace_back("ANY");
//                for (auto &np : nodeTypes) {
//                    nodeTypePattern.push_back(np);
//                    break;
//                }
            } else {
                edgeTypePattern.push_back(pathIns[i]);
            }
        }
        if (abandon)
            continue;

        ///Expand the node
        size_t edge_start = gm.nodes[current_node];
        size_t adj_sz = gm.nodes[current_node + 1] - edge_start;

        bool handledHD = hdNid2EffectNeighbors.count(current_node);
        for (int i = 0; i < adj_sz; ++i) {
            uint64_t edge;
            if (handledHD) {
                //avoid scan the neighbors of hd nodes again and again
                if (i >= hdNid2EffectNeighbors[current_node].size()) break;
                edge = hdNid2EffectNeighbors[current_node][i];
            } else {
                edge = gm.edges[edge_start + i];
            }

            uint32_t eid = extractHigh32bits(edge); // edge
            uint32_t end_node = extractLow32bits(edge); // end node
            double sz_exp = nodeTypePattern.size();

            if (end_node == target) continue; //inside node cannot be either contextNode or target node
            if (current_node == end_node) continue; //no self-cycle

            //commented for top-k search
//            if (!nodeScore.count(end_node) || nodeScore[end_node] < relevancy_limit) {
//                continue; // the node score is too low. cancel path instance level relevancy limit
//            }

            if (!nodeScore.count(end_node)) {
                continue; // the node score is too low. cancel path instance level relevancy limit
            } else if (nodeScore[end_node] < current_kth_relevance_score()) {
                continue;
            }

            if (adj_sz >= degree_threshold && !handledHD) {
                hdNid2EffectNeighbors[current_node].push_back(edge);
            }

            double end_node_score = nodeScore[end_node] * pow(base_e, -sz_exp);
            double path_score = end_node_score < curr_pair.first ? end_node_score : curr_pair.first;

            //commented for top-k search
//            if (path_score < relevancy_limit) {
//                continue; // cancel path instance level relevancy
//            }
            if (path_score < current_kth_relevance_score())
                continue;

            //check length.

            ///Check possible repetitions
            vector<uint32_t> tmpEdgeTypes = edgeTypePattern;
            tmpEdgeTypes.push_back(eid);

            auto sigRep = pathPattern2SigWNodes(nodeTypePattern, tmpEdgeTypes);
            if (endNode2pathSig.count(end_node) && endNode2pathSig[end_node].count(sigRep))
                continue;
            endNode2pathSig[end_node].insert(sigRep);

            /// use reverse info
            test1 = getTime();
            backWardFromNode_DFS_Reverse(gm, tmpEdgeTypes, visitedEdgeSeq,
                    sz_exp, end_node, target, hop_num, path_score);
            // backWardFromNode_Naive(gm, tmpEdgeTypes, visitedEdgeSeq, sz_exp, end_node, target, hop_num);
            test2 = getTime();
            backwardFromNodeTime += getInterval(test1, test2);
            if (topk_relevant_pq.size() < top_pattern_k)
                backwardTimeBeforetopkCollected += getInterval(test1, test2);

            // if (backwardFromNodeTime > timeLimit * 1000) return;
            if (contextNodes.count(end_node)) continue;

            //with index, we only need one condition


            //below for handling low degree node
            if (sz_exp + 1 + 2 > path_length_limit) {
                continue;
            }

            double possible_largest =
                    smallest_one_hop * pow(base_e, -(sz_exp + 1)); // no high degree node Or with index
            possible_largest = path_score < possible_largest ? path_score : possible_largest;

            //commented for top-k search
//            if (possible_largest < relevancy_limit) {
//                continue; //cancel path instance level relevancy limit
//            }
            if (possible_largest < current_kth_relevance_score()) {
                continue;
            }

            vector<uint32_t> tmpIns = pathIns;
            tmpIns.push_back(eid);
            tmpIns.push_back(end_node);
            pathQueue.push(make_pair(path_score, make_pair(possible_largest, gm.pathPattern2sig(tmpIns))));
        }
//        auto inloop_e_time = getTime();
//        printf("Inloop time = %.3lf ms, nid = %d, Adj_sz = %ld, iri = %s\n", getInterval(inloop_time, inloop_e_time),
//               current_node, adj_sz, gm.nid2IRI[current_node].c_str());
    }


    uint64_t inner_e = getTime();
    InnerProcessingTime = getInterval(inner_s, inner_e);
    ///print final info
//    cout << "#### Output global top-k exceptional facts" << endl;
//    if (res_pq.empty() || res_pq.top().excepScore < 0.7)
//        return;
    printResQueue(gm, target, res_pq);
//    cout << "\n#### Output global Diverse top-k exceptional facts" << endl;
//    printResQueue(gm, target, diverseTopk);
//
//    printf("Inner time = %.3lf ms.\n", getInterval(inner_s, inner_e));

}


void ExceptionalFactMining::backWardFromNode_DFS_Reverse(GraphManager &gm, std::vector<uint32_t> &tmpEdgeTypes,
                                                         std::unordered_set<std::string> &visitedEdgeSeq, double sz_exp,
                                                         int node2add, int target, int hop_num,
                                                         double path_score) {
    for (int j = 0; j < hop_num; ++j) {
        if ((int) sz_exp + 1 + (1 + j) > path_length_limit) break;
        ///check for low degree nodes
        if (!hopNid2Parent_vec.empty() && hopNid2Parent_vec[j].count(node2add)) {
            for (const auto &s: hopNid2Parent_vec[j][node2add]) {
                //judge score first
                double s_score = s.first * pow(base_e, -sz_exp);

                // commented for using top-k relevance
//                if (s_score < relevancy_limit) {
//                    //opt_1++;
//                    /// Important!!! must use "break" here
//                    /// because there may exist other nodes with the same edge seq but has higher score
//                    break; // path instance level relevancy
//                }

                if (s_score < current_kth_relevance_score()) {
                    break;
                }


                vector<uint32_t> reverseEdgePattern;
                vector<unordered_set<uint32_t>> levelFrontiers;

                auto contxtNid = getReversePathType(gm, tmpEdgeTypes, s.second, reverseEdgePattern);
                auto edgeSeqSig = reverseEdgeSeqSig(gm, reverseEdgePattern);

                if (visitedEdgeSeq.count(edgeSeqSig)) {
                    continue;
                }
                visitedEdgeSeq.insert(edgeSeqSig);
                levelFrontiers.resize(reverseEdgePattern.size() + 1);
                levelFrontiers.back().insert(contxtNid);

                nodeSeqLabels_DFS.clear();
                double pScore = pathPatternFromContextNode_DFS(gm, reverseEdgePattern,
                                                               levelFrontiers,
                                                               (int) reverseEdgePattern.size(),
                                                               target, contxtNid);

                pScore = pScore < path_score ? pScore : path_score;

                if (pScore > 0 && levelFrontiers.front().size() >= num_peerEntity) {
                    auto t1 = getTime();
                    for (int i = 1; i < levelFrontiers.size() - 1; ++i) {
                        if (levelFrontiers[i].size() == 1) {
                            nodeSeqLabels_DFS[i - 1].append("(");
                            for (auto &f : levelFrontiers[i])
                                nodeSeqLabels_DFS[i - 1].append(gm.nid2IRI[f]);
                            nodeSeqLabels_DFS[i - 1].append(")");
                        }
                    }

                    emitResult(gm, edgeSeqSig, levelFrontiers.front(), target, contxtNid, pScore);
                    //emitResultJaccardSim(gm, edgeSeqSig, levelFrontiers.front(), target, contxtNid, pScore);
                    auto t2 = getTime();
                    emitResTime += getInterval(t1, t2);
                }

            }
        }
    }
}


void ExceptionalFactMining::oneHopNeighborFind(GraphManager &gm, uint32_t target) {
    unordered_set<uint32_t> visitedEid;
    int edge_start = gm.nodes[target];
    int adj_sz = gm.nodes[target + 1] - edge_start;

    for (int i = 0; i < adj_sz; ++i) {
        auto edge = gm.edges[edge_start + i];
        auto end_node = extractLow32bits(edge);
        auto eid = extractHigh32bits(edge);
        if (visitedEid.count(eid)) continue;
        visitedEid.insert(eid);

        auto reverseEid = gm.reverseEdgeTypeId(eid);
        unordered_set<uint32_t> _fset;
        int _ePos = gm.snidEid2pos[combine2u32(end_node, reverseEid)];
        if (_ePos == -1) {
            cout << "Impoossible!" << endl;
            continue;
        }

        for (int j = _ePos; j < gm.nodes[end_node + 1]; ++j) {
            auto reid = extractHigh32bits(gm.edges[j]);
            if (reid != reverseEid) break;

            _fset.insert(extractLow32bits(gm.edges[j]));

        }

        if (_fset.size() < num_peerEntity) continue;
        nodeSeqLabels_DFS.clear();
        nodeSeqLabels_DFS.push_back(gm.nid2IRI[end_node]);
        string edgeSeqSig = to_string(eid);
        emitResult(gm, edgeSeqSig, _fset, target, 0, 1);

    }

    cout << "#### Output global top-k exceptional facts" << endl;
    printResQueue(gm, target, res_pq);
    cout << "\n#### Output global Diverse top-k exceptional facts" << endl;
    printResQueue(gm, target, diverseTopk);
}


void ExceptionalFactMining::backWardFromNode_Blind_EndingNode(GraphManager &gm, std::vector<uint32_t> &tmpEdgeTypes,
                                                   std::unordered_set<std::string> &visitedEdgeSeq, int node2add,
                                                   int target) {
    if (node2add == target) return;
    vector<unordered_set<uint32_t>> levelFrontiers;
    auto reverseEdgePattern = gm.reversePatternSequence(tmpEdgeTypes);
    auto edgeSeqSig = reverseEdgeSeqSig(gm, reverseEdgePattern);
    if (visitedEdgeSeq.count(edgeSeqSig)) {
        return;
    }
    visitedEdgeSeq.insert(edgeSeqSig);
    levelFrontiers.resize(reverseEdgePattern.size() + 1);
    levelFrontiers.back().insert(node2add);

    //BFS search for peer entities
    nodeSeqLabels_DFS.clear();
    nodeSeqLabels_DFS.resize(reverseEdgePattern.size());
    for (int i = 0; i < reverseEdgePattern.size(); ++i) {
        auto eid = reverseEdgePattern[i];
        int levelF_idx = (int) levelFrontiers.size() - i - 1 - 1;
        bool setLabel = false;
        for (auto f : levelFrontiers[levelF_idx + 1]) {
            auto _key = combine2u32(f, eid);
            if (!gm.snidEid2pos.count(_key)) continue;
            if (!setLabel) {
                string _l;
                gm.getOneLabelRandom(f, _l);
                nodeSeqLabels_DFS[levelF_idx] = _l;
                setLabel = true;
            }

            auto edgePos = gm.snidEid2pos[_key];
            for (int j = edgePos; j < gm.nodes[f + 1]; ++j) {
                auto _node = extractLow32bits(gm.edges[j]);
                levelFrontiers[levelF_idx].insert(_node);
            }
        }
    }
    nodeSeqLabels_DFS.back() = gm.nid2IRI[node2add];

    if (levelFrontiers.front().size() >= num_peerEntity) {
        auto t1 = getTime();
        double pScore = 1;
        emitResult_Blind(gm, edgeSeqSig, levelFrontiers.front(), target, node2add, pScore);
        auto t2 = getTime();
        emitResTime += getInterval(t1, t2);
    }

}


void ExceptionalFactMining::backWardFromNode_Blind_EndingType(GraphManager &gm, std::vector<uint32_t> &tmpEdgeTypes,
                                                              std::unordered_set<std::string> &visitedEdgeSeq, int node2add,
                                                              int target) {
    if (node2add == target) return;
    vector<unordered_set<uint32_t>> levelFrontiers;
    auto reverseEdgePattern = gm.reversePatternSequence(tmpEdgeTypes);
    auto edgeSeqSig = reverseEdgeSeqSig(gm, reverseEdgePattern);
    if (visitedEdgeSeq.count(edgeSeqSig)) {
        return;
    }
    visitedEdgeSeq.insert(edgeSeqSig);
    levelFrontiers.resize(reverseEdgePattern.size() + 1);

    string node2addType;
    gm.getOneLabelRandom(node2add, node2addType);
    if (node2addType == "ANY") {
        levelFrontiers.back() = gm.typeId2Count[reverseEdgePattern.front()];
    } else {
        for (auto nf : gm.type2nid[node2addType])
            if (gm.typeId2Count[reverseEdgePattern.front()].count(nf))
                levelFrontiers.back().insert(nf);
    }

    //BFS search for peer entities
    nodeSeqLabels_DFS.clear();
    nodeSeqLabels_DFS.resize(reverseEdgePattern.size());
    for (int i = 0; i < reverseEdgePattern.size(); ++i) {
        auto eid = reverseEdgePattern[i];
        int levelF_idx = (int) levelFrontiers.size() - i - 1 - 1;
        bool setLabel = false;
        for (auto f : levelFrontiers[levelF_idx + 1]) {
            auto _key = combine2u32(f, eid);
            if (!gm.snidEid2pos.count(_key)) continue;
            if (!setLabel) {
                string _l;
                gm.getOneLabelRandom(f, _l);
                nodeSeqLabels_DFS[levelF_idx] = _l;
                setLabel = true;
            }

            auto edgePos = gm.snidEid2pos[_key];
            for (int j = edgePos; j < gm.nodes[f + 1]; ++j) {
                auto _node = extractLow32bits(gm.edges[j]);
                levelFrontiers[levelF_idx].insert(_node);
            }

        }
    }

    if (levelFrontiers.front().size() >= num_peerEntity) {
        auto t1 = getTime();
        double pScore = 1;
        emitResult_Blind(gm, edgeSeqSig, levelFrontiers.front(), target, node2add, pScore);
        auto t2 = getTime();
        emitResTime += getInterval(t1, t2);
    }

}

void ExceptionalFactMining::enumPathPatternBlind(GraphManager &gm, uint32_t target) {
    uint64_t inner_s = getTime();
    //khopReverseExpandingFromContext(gm, target, hop_num);

    unordered_map<size_t, unordered_map<uint32_t, string>> peerSz2firstPthIns;

    /// The string refers to the path instance.
    /// From path instance, we can infer the corresponding path pattern.
    /// The second element is the signature of path instance.
    typedef pair<double, pair<double, string>> patternSeqScore;
    /// By default, this is a max heap that is needed.
    priority_queue<patternSeqScore, vector<patternSeqScore>> pathQueue;
    vector<uint32_t> initPath;
    initPath.push_back(target);
    string _sig = gm.pathPattern2sig(initPath);
    double _initScore = 1;
    pathQueue.push(make_pair(_initScore, make_pair(1, _sig)));

    /// For repetition checking
    unordered_map<uint32_t, unordered_set<string>> endNode2pathSig;
    unordered_set<string> visitedEdgeSeq;

    while (!pathQueue.empty()) {
        auto tLimit = getTime();
        if (getInterval(inner_s, tLimit) >= TimeLimit) {
            InnerProcessingTime = TimeLimit;
            break; // timeout still needs to out put answer
        }

        auto curr_pair = pathQueue.top();
        pathQueue.pop();
        auto pathInsSig = curr_pair.second;

        ///Get path patterns
        auto pathIns = gm.sig2pathPattern(pathInsSig.second);
        vector<string> nodeTypePattern;
        vector<uint32_t> edgeTypePattern;
        uint32_t current_node = pathIns.back();

        for (int i = 0; i < pathIns.size(); ++i) {
            if (i % 2 == 0) {
                //it is a node
                auto nodeTypes = gm.nid2types[pathIns[i]]; // node type
                // if this is the case, then the node has no type. We fix this node
                if (nodeTypes.empty())
                    nodeTypePattern.emplace_back("ANY");
                for (auto &np : nodeTypes) {
                    nodeTypePattern.push_back(np);
                    break;
                }
            } else {
                edgeTypePattern.push_back(pathIns[i]);
            }
        }

        double sz_exp = nodeTypePattern.size();
        if (sz_exp >= path_length_limit) continue;

        ///Expand the node
        size_t edge_start = gm.nodes[current_node];
        size_t adj_sz = gm.nodes[current_node + 1] - edge_start;

        for (int i = 0; i < adj_sz; ++i) {
            auto edge = gm.edges[edge_start + i];
            uint32_t eid = extractHigh32bits(edge);
            uint32_t end_node = extractLow32bits(edge);

            if (end_node == target) continue; //inside node cannot be either contextNode or target node
            if (current_node == end_node) continue; //no self-cycle

            ///Check possible repetitions
            vector<uint32_t> tmpEdgeTypes = edgeTypePattern;
            tmpEdgeTypes.push_back(eid);

            auto sigRep = pathPattern2SigWNodes(nodeTypePattern, tmpEdgeTypes);


            if (endNode2pathSig.count(end_node) && endNode2pathSig[end_node].count(sigRep))
                continue;
            endNode2pathSig[end_node].insert(sigRep);


            auto test1 = getTime();
            backWardFromNode_Blind_EndingNode(gm, tmpEdgeTypes, visitedEdgeSeq, end_node, target);
            backWardFromNode_Blind_EndingType(gm, tmpEdgeTypes, visitedEdgeSeq, end_node, target);
            auto test2 = getTime();
            backwardFromNodeTime += getInterval(test1, test2);


            vector<uint32_t> tmpIns = pathIns;
            tmpIns.push_back(eid);
            tmpIns.push_back(end_node);
            pathQueue.push(make_pair(1.0 / (double) tmpIns.size(),
                                     make_pair(1.0 / (double) tmpIns.size(), gm.pathPattern2sig(tmpIns))));
        }
    }


    uint64_t inner_e = getTime();

    ///print final info
    cout << "\n#### Output global top-k exceptional facts" << endl;
    printResQueue(gm, target, res_pq);
//    cout << "\n#### Output global Diverse top-k exceptional facts" << endl;
//    printResQueue(gm, target, diverseTopk);

    auto tmp_InnerProcessingTime = getInterval(inner_s, inner_e);
    InnerProcessingTime = InnerProcessingTime < tmp_InnerProcessingTime? InnerProcessingTime : tmp_InnerProcessingTime;
}





//void
//ExceptionalFactMining::enumPathPatternWithReverseExpansion_DFS_Reverse_MemorizeHD(GraphManager &gm,
//                                                                                  uint32_t target,
//                                                                                  int hop_num) {
//
//    uint64_t inner_s = getTime();
//    khopReverseExpandingFromContext_StrOnly(gm, target, hop_num);
//
//    unordered_map<size_t, unordered_map<uint32_t, string>> peerSz2firstPthIns;
//
//    /// The string refers to the path instance.
//    /// From path instance, we can infer the corresponding path pattern.
//    /// The second element is the signature of path instance.
//    typedef pair<double, pair<double, string>> patternSeqScore;
//    /// By default, this is a max heap that is needed.
//    priority_queue<patternSeqScore, vector<patternSeqScore>> pathQueue;
//    vector<uint32_t> initPath;
//    initPath.push_back(target);
//    string _sig = gm.pathPattern2sig(initPath);
//    double _initScore = 1;
//    pathQueue.push(make_pair(_initScore, make_pair(1, _sig)));
//
//    /// For repetition checking
//    unordered_map<uint32_t, unordered_set<string>> endNode2pathSig;
//    unordered_set<string> visitedEdgeSeq;
//
//    /// For preventing the path patterns leading to unchanged peer entities
//    unordered_map<uint32_t, unordered_map<string, unordered_set<uint32_t>>> node2pattern2peerEntities;
//
//
//    /// need to handle the single target node without expansion
//    vector<uint32_t> emptyEdgeTypes;
//    auto test1 = getTime();
//    backWardFromNode_DFS_Reverse_MemorizeHD(gm, emptyEdgeTypes, visitedEdgeSeq, 0, target, target, hop_num);
////    backWardFromNode_Naive(gm, emptyEdgeTypes, visitedEdgeSeq, 0, target, target, hop_num);
//    auto test2 = getTime();
//    backwardFromNodeTime += getInterval(test1, test2);
//
//    while (!pathQueue.empty()) {
//        auto curr_pair = pathQueue.top();
//        pathQueue.pop();
//        auto pathInsSig = curr_pair.second;
//
//        ///Get path patterns
//        auto pathIns = gm.sig2pathPattern(pathInsSig.second);
//        vector<string> nodeTypePattern;
//        vector<uint32_t> edgeTypePattern;
//        uint32_t current_node = pathIns.back();
//
//        for (int i = 0; i < pathIns.size(); ++i) {
//            if (i % 2 == 0) {
//                //it is a node
//                auto nodeTypes = gm.nid2types[pathIns[i]]; // node type
//                // if this is the case, then the node has no type. We fix this node
//                if (nodeTypes.empty())
//                    nodeTypePattern.emplace_back("ANY");
//                for (auto &np : nodeTypes) {
//                    nodeTypePattern.push_back(np);
//                    break;
//                }
//            } else {
//                edgeTypePattern.push_back(pathIns[i]);
//            }
//        }
//
//        ///Expand the node
//        size_t edge_start = gm.nodes[current_node];
//        size_t adj_sz = gm.nodes[current_node + 1] - edge_start;
//
//        bool handledHD = hdNid2EffectNeighbors.count(current_node);
//        for (int i = 0; i < adj_sz; ++i) {
//            uint64_t edge;
//            if (handledHD) {
//                //avoid scan the neighbors of hd nodes again and again
//                if (i >= hdNid2EffectNeighbors[current_node].size()) break;
//                edge = hdNid2EffectNeighbors[current_node][i];
//            } else {
//                edge = gm.edges[edge_start + i];
//            }
//
//
//            uint32_t eid = extractHigh32bits(edge);
//            uint32_t end_node = extractLow32bits(edge);
//            double sz_exp = nodeTypePattern.size();
//
//            if (end_node == target) continue; //inside node cannot be either contextNode or target node
//            if (current_node == end_node) continue; //no self-cycle
//            if (!nodeScore.count(end_node) || nodeScore[end_node] < relevancy_limit) {
//                continue; // the node score is too low. cancel path instance level relevancy limit
//            }
//            if (adj_sz >= degree_threshold && !handledHD) {
//                hdNid2EffectNeighbors[current_node].push_back(edge);
//            }
//
//
//            double end_node_score = nodeScore[end_node] * pow(base_e, -sz_exp);
//            double path_score = end_node_score < curr_pair.first ? end_node_score : curr_pair.first;
//            if (path_score < relevancy_limit) {
//                continue; // cancel path instance level relevancy
//            }
//            //check length.
//
//            ///Check possible repetitions
//            vector<uint32_t> tmpEdgeTypes = edgeTypePattern;
//            tmpEdgeTypes.push_back(eid);
//
//            auto sigRep = pathPattern2SigWNodes(nodeTypePattern, tmpEdgeTypes);
//            if (endNode2pathSig.count(end_node) && endNode2pathSig[end_node].count(sigRep))
//                continue;
//            endNode2pathSig[end_node].insert(sigRep);
//
//            /// use reverse info
//            test1 = getTime();
//            backWardFromNode_DFS_Reverse_MemorizeHD(gm, tmpEdgeTypes, visitedEdgeSeq, sz_exp, end_node, target,
//                                                    hop_num);
//            // backWardFromNode_Naive(gm, tmpEdgeTypes, visitedEdgeSeq, sz_exp, end_node, target, hop_num);
//            test2 = getTime();
//            backwardFromNodeTime += getInterval(test1, test2);
//
//            if (contextNodes.count(end_node)) continue;
//
//            //with index, we only need one condition
//
//
//            //below for handling low degree node
//            if (sz_exp + 1 + 2 > path_length_limit) {
//                continue;
//            }
//
//            double possible_largest =
//                    smallest_one_hop * pow(base_e, -(sz_exp + 1)); // no high degree node Or with index
//            possible_largest = path_score < possible_largest ? path_score : possible_largest;
//            if (possible_largest < relevancy_limit) {
//                continue; //cancel path instance level relevancy limit
//            }
//
//            vector<uint32_t> tmpIns = pathIns;
//            tmpIns.push_back(eid);
//            tmpIns.push_back(end_node);
//            pathQueue.push(make_pair(path_score, make_pair(possible_largest, gm.pathPattern2sig(tmpIns))));
//        }
//    }
//
//
//    uint64_t inner_e = getTime();
//
//    ///print final info
//    cout << "\n#### Output global top-k exceptional facts" << endl;
//    printResQueue(gm, res_pq, node2pattern2peerEntities);
//    cout << "\n#### Output global Diverse top-k exceptional facts" << endl;
//    printResQueue(gm, diverseTopk, node2pattern2peerEntities);
//
//    printf("Inner time = %.3lf ms.\n", getInterval(inner_s, inner_e));
//
//}
//
//void ExceptionalFactMining::backWardFromNode_DFS_Reverse_MemorizeHD(GraphManager &gm,
//                                                                    std::vector<uint32_t> &tmpEdgeTypes,
//                                                                    std::unordered_set<std::string> &visitedEdgeSeq,
//                                                                    double sz_exp, int node2add, int target,
//                                                                    int hop_num) {
//
//    for (int j = 0; j < hop_num; ++j) {
//        if ((int) sz_exp + 1 + (1 + j) > path_length_limit) break;
//        ///check for low degree nodes
//        if (!hopNid2Parent_vec.empty() && hopNid2Parent_vec[j].count(node2add)) {
//            for (const auto &s: hopNid2Parent_vec[j][node2add]) {
//                //judge score first
//                double s_score = s.first * pow(base_e, -sz_exp);
//                if (s_score < relevancy_limit) {
//                    //opt_1++;
//                    /// Important!!! must use "break" here
//                    /// because there may exist other nodes with the same edge seq but has higher score
//                    break; // path instance level relevancy
//                }
//
////                if (s.second == "82 83 3198") {
////                    cout << "sss1111sss" << endl;
////                }
////                if (s.second == "54 55 3198") {
////                    cout << "ssssss" << endl;
////                }
//                vector<uint32_t> reverseEdgePattern;
//                vector<unordered_set<uint32_t>> levelFrontiers;
//
//                auto contxtNid = getReversePathType(gm, tmpEdgeTypes, s.second, reverseEdgePattern);
//                auto edgeSeqSig = reverseEdgeSeqSig(gm, reverseEdgePattern);
//
//                if (visitedEdgeSeq.count(edgeSeqSig)) {
//                    continue;
//                }
//                visitedEdgeSeq.insert(edgeSeqSig);
//
//                levelFrontiers.resize(reverseEdgePattern.size() + 1);
//                levelFrontiers.back().insert(contxtNid);
////                for (int k = tmpEdgeTypes.size(); k < levelFrontiers.size() - 1; ++k) {
////                    int k_idx = (int) levelFrontiers.size() - k - 2;
////                    levelFrontiers[k] = pattern2rightFrontiers_vec[s.second][k_idx];
////                }
//                auto time1 = getTime();
//
//                nodeSeqLabels_DFS.clear();
//                double pScore = pathPatternFromContextNode_DFS_MemorizeHD(gm, reverseEdgePattern,
//                                                                          levelFrontiers,
//                                                                          (int) reverseEdgePattern.size(),
//                                                                          target, contxtNid);
//
//                auto time2 = getTime();
//                pathPatternFromContextNodeTime += getInterval(time1, time2);
//                if (pScore > 0 && levelFrontiers.front().size() >= num_peerEntity) {
//                    emitResult(gm, edgeSeqSig, levelFrontiers.front(), target, contxtNid, pScore);
//                }
//
//            }
//        }
//    }
//}

//bool ExceptionalFactMining::DFS_pattern_score_MemorizeHD(GraphManager &gm,
//                                                         std::vector<uint32_t> &reversedEdgeSeq,
//                                                         int edge_ptr,
//                                                         std::vector<std::unordered_set<uint32_t>> &levelFrontiers,
//                                                         std::vector<std::unordered_set<uint32_t>> &levelHDNodes,
//                                                         std::vector<double> &levelScore,
//                                                         std::vector<double> &numFrontier,
//                                                         uint32_t cf, uint32_t source, bool &discard,
//                                                         bool parentHD, uint32_t parentHDId) {
//    if (discard) return false;
//    if (edge_ptr == reversedEdgeSeq.size()) {
//        levelFrontiers[0].insert(cf); // last edge
//        if (!contextNodes.count(cf)) {
//            numFrontier[0]++;
//            if (nodeScore.count(cf))
//                levelScore[0] += nodeScore[cf];
//        }
//        //test limit num
//        if (!naiveDFS)
//            if (testFrontierSzLimit(gm, (int) numFrontier[0], reversedEdgeSeq.back())) {
//                discard = true;
//                return false;
//            }
//        return true;
//    }
//    int levelFrontier_idx = (int) reversedEdgeSeq.size() - edge_ptr - 1;
//    int edge_pos = gm.firstEdgePos(cf, reversedEdgeSeq[edge_ptr]);
//    if (edge_pos == -1) return false;
//    bool inserted = false;
//    for (int i = edge_pos; i < gm.nodes[cf + 1]; ++i) {
//        auto eid = extractHigh32bits(gm.edges[i]);
//        if (eid != reversedEdgeSeq[edge_ptr]) break;
//        uint32_t end_node = extractLow32bits(gm.edges[i]);
//        //one edge look forward
//        if (!naiveDFS)
//            if (edge_ptr < reversedEdgeSeq.size() - 1) {
//                uint32_t oneMoreEdge = reversedEdgeSeq[edge_ptr + 1];
//                if (!gm.typeId2Count[oneMoreEdge].count(end_node)) {
//                    continue;
//                }
//            }
//
//
//        bool valid = levelFrontiers[levelFrontier_idx].count(end_node);
//        if (!valid)
//            valid = DFS_pattern_score_MemorizeHD(gm, reversedEdgeSeq, edge_ptr + 1,
//                                                 levelFrontiers, levelHDNodes, levelScore, numFrontier, end_node,
//                                                 source, discard);
//        if (discard) {
//            return false;
//        }
//        if (valid) {
//            if (!inserted) {
//                levelFrontiers[levelFrontier_idx + 1].insert(cf);
//                inserted = true;
//                if (!contextNodes.count(cf)) {
//                    numFrontier[levelFrontier_idx + 1]++;
//                    if (nodeScore.count(cf))
//                        levelScore[levelFrontier_idx + 1] += nodeScore[cf];
//
//                    //later for printing labels
//                    if (edge_ptr != 0 && edge_ptr + 1 + nodeSeqLabels_DFS.size() == reversedEdgeSeq.size()) {
//                        string _label;
//                        gm.getOneLabelRandom(cf, _label);
//                        nodeSeqLabels_DFS.push_back(_label);
//                    }
//
//                }
//                if (!naiveDFS)
//                    if (edge_ptr > 0)
//                        if (testFrontierSzLimit(gm, (int) numFrontier[levelFrontier_idx + 1],
//                                                reversedEdgeSeq[edge_ptr - 1])) {
//                            discard = true;
//                            return false;
//                        }
//
//            }
//        }
//    }
//    return inserted;
//}


//double ExceptionalFactMining::pathPatternFromContextNode_DFS_MemorizeHD(GraphManager &gm,
//                                                                        std::vector<uint32_t> &reversedEdgeSeq,
//                                                                        std::vector<std::unordered_set<uint32_t>> &levelFrontiers,
//                                                                        int mid_idx, int target, int source) {
//    vector<double> levelScore(reversedEdgeSeq.size() + 1);
//    vector<double> numFrontier(reversedEdgeSeq.size() + 1);
//    vector<unordered_set<uint32_t>> levelHDNodes(reversedEdgeSeq.size() + 1); // actually. No need to plus 1 to the size
//
//    levelScore.back() += nodeScore[source];
//    numFrontier.back()++;
//
//    bool discard = false;
//    int rest_sz = (int) reversedEdgeSeq.size() - mid_idx;
//    int levelF_idx = mid_idx - 1;
//    for (auto cf : levelFrontiers[levelF_idx + 1]) {
//        bool valid = DFS_pattern_score_MemorizeHD(gm, reversedEdgeSeq, rest_sz,
//                                                  levelFrontiers, levelHDNodes, levelScore, numFrontier,
//                                                  cf, source,
//                                                  discard,
//                                                  false, cf); // here the context node we do not care
//        if (discard) {
//            opt_1++;
//            return -1;
//        }
//        if (valid) {
//            levelFrontiers[levelF_idx + 1].insert(cf);
//        }
//    }
//
//
//    if (levelFrontiers.front().size() < num_peerEntity)
//        return -1;
//    double s = 1;
//    for (int i = 0; i < levelScore.size() - 1; ++i) {
//        if (numFrontier[i] == 0) return -1;
//        double tmp_s = levelScore[i] / numFrontier[i];
//        tmp_s *= pow(base_e, -i);
//        s = s < tmp_s ? s : tmp_s;
//        if (s < relevancy_limit) return -1;
//    }
//    nodeSeqLabels_DFS.push_back(gm.nid2IRI[source]);
//    return s;
//}
//
