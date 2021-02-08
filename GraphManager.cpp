//
// Created by yangyueji on 7/15/19.
//

#include "GraphManager.h"
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <queue>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/serialization/vector.hpp>

using namespace std;

int GraphManager::firstEdgePos(uint32_t start_node, uint32_t eid) const {
    uint64_t edgeValue = combine2u32(eid, 0);
    int _min = nodes[start_node];
    int _max = nodes[start_node + 1];
    while (_min < _max) {
        int _d = (_min + _max) / 2;
        if (edges[_d] > edgeValue) {
            _max = _d;
        } else {
            _min = _d + 1;
        }
    }
    if (extractHigh32bits(edges[_min]) == eid)
        return _min;
    return -1;

}


void GraphManager::readEdges(const char *filename) {

    unordered_map<string, int> etypes2id;
    unordered_map<string, uint32_t> IRI2nid;

    vector<pair<uint64_t, string>> tmp_edges;
    ifstream ifs(filename);
    if (!ifs.is_open()) {
        cout << "Invalid filepath\n";
        return;
    }

    string line;
    bool v1 = false, v2 = false;
    while (getline(ifs, line)) {
        istringstream iss(line);
        uint32_t sid, eid;
        string token;
        getline(iss, token, ',');
        if (IRI2nid.count(token))
            sid = IRI2nid[token];
        else {
            sid = (uint32_t) IRI2nid.size();
            IRI2nid[token] = sid;
        }

        if (!v1 && token == "http://www.wikidata.org/entity/Q22686") {
            cout << "trump, " << sid << endl;
            v1 = true;
        }

        if (!v2 && token == "http://www.wikidata.org/entity/Q6294") {
            v2 = true;
            cout << "Clinton, " << sid << endl;
        }

        getline(iss, token, ',');
        if (IRI2nid.count(token))
            eid = IRI2nid[token];
        else {
            eid = (uint32_t) IRI2nid.size();
            IRI2nid[token] = eid;
        }

        if (sid == eid) continue;

        getline(iss, token, ',');
        string reversedType = "r-" + token;
        if (!etypes2id.count(token)) {
            int sz = etypes2id.size();
            etypes2id[token] = sz;
            etypes2id[reversedType] = sz + 1;
        }
        tmp_edges.emplace_back(make_pair(combine2u32(sid, eid), token));
        tmp_edges.emplace_back(make_pair(combine2u32(eid, sid), reversedType));
    }

    ifs.close();

    sort(tmp_edges.begin(), tmp_edges.end());
    //remove repeatitions
    auto iter = std::unique(tmp_edges.begin(), tmp_edges.end());
//    cout << extractHigh32bits(tmp_edges[303010].first) << "," << extractLow32bits(tmp_edges[303010].first) << ","
//         << tmp_edges[303010].second << endl;
//    cout << extractHigh32bits(tmp_edges[303011].first) << "," << extractLow32bits(tmp_edges[303011].first) << ","
//         << tmp_edges[303011].second << endl;
    tmp_edges.resize(std::distance(tmp_edges.begin(), iter));

    node_sz = IRI2nid.size();
    edge_sz = tmp_edges.size();
    nodes = new uint32_t[node_sz + 1];
    edges = new uint64_t[edge_sz];
    edgesSpreadRatio = new double[edge_sz];
//    uint32_t last_index = 0;
    uint32_t iter_j = 0;
    for (int i = 0; i < node_sz; ++i) {
        uint32_t startNode = extractHigh32bits(tmp_edges[iter_j].first);
        if (i <= startNode) {
            nodes[i] = iter_j;
            continue;
        }

        if (iter_j == edge_sz) {
            nodes[i] = (uint32_t) edge_sz;
            continue;
        }
        while (iter_j < edge_sz && extractHigh32bits(tmp_edges[iter_j].first) ==
                                   startNode) {
            iter_j++;
        }
        nodes[i] = iter_j;
//        last_index = iter_j;
    }
    nodes[node_sz] = (uint32_t) edge_sz;

//    cout << etypes2id["diplomatic relation"] << endl;
    /**Set edges type*/
    for (int i = 0; i < edge_sz; ++i) {
        int tid = etypes2id[tmp_edges[i].second];
        uint32_t eid = extractLow32bits(tmp_edges[i].first);
        edges[i] = combine2u32((uint32_t) tid, eid);
    }

    typeId2Count.resize(etypes2id.size());
    for (uint32_t i = 0; i < node_sz; ++i) {
        //re-sort neighbors by edge types;
        int adj_sz = nodes[i + 1] - nodes[i];
        if (adj_sz == 0) continue;

        sort(edges + nodes[i], edges + nodes[i] + adj_sz);
        unordered_map<uint32_t, double> tid2count;
        for (int j = 0; j < adj_sz; ++j) {
            auto e = edges[nodes[i] + j];
            uint32_t tid = extractHigh32bits(e);
//            if (i == 3309) {
//                cout << tid << "," << extractLow32bits(e) << endl;
//            }
            typeId2Count[tid].insert(i);
            if (tid2count.count(tid))
                tid2count[tid] += 1;
            else
                tid2count[tid] = 1;
        }

        for (int j = 0; j < adj_sz; ++j) {
            auto e = edges[nodes[i] + j];
            uint32_t tid = extractHigh32bits(e);
            edgesSpreadRatio[nodes[i] + j] = 1.0 / tid2count[tid];
        }


    }

    typeId2Name.resize(etypes2id.size());
    for (auto &kv_pair : etypes2id)
        typeId2Name[kv_pair.second] = kv_pair.first;

    nid2IRI.resize(IRI2nid.size());
    for (auto &kv_pair : IRI2nid)
        nid2IRI[kv_pair.second] = kv_pair.first;



    printf("reading finished!\n");
}

void GraphManager::readNodeTypes(const char *filename) {
    ifstream ifs(filename);
    if (!ifs.is_open()) {
        cout << "Invalid File Path" << endl;
    }

    nid2types.resize(nid2IRI.size());
    unordered_map<string, uint32_t> iri2node;
    for (uint32_t nid = 0; nid < nid2IRI.size(); ++nid) {
        iri2node[nid2IRI[nid]] = nid;
    }


    string line;
    while (getline(ifs, line)) {
        istringstream iss(line);
        string iri, t;
        unordered_set<string> types;
        getline(iss, iri, ','); // node IRI
        if (!iri2node.count(iri)) continue; // node without edges are removed
        while (getline(iss, t, ',')) {
            if (t.empty()) continue;
            types.insert(t);
            type2nid[t].insert(iri2node[iri]);
        }
        nid2types[iri2node[iri]] = types;

    }

    ifs.close();
}

void GraphManager::resetNodeTypes() {
    printf("Get most specific node type as starting point of path patterns.\n");
//    nid2types.clear();
//    type2nid.clear();
    auto s1 = getTime();
    unordered_map<string, uint32_t> type2typeId;
    for (uint32_t i = 0; i < typeId2Name.size(); ++i) {
        type2typeId[typeId2Name[i]] = i;
    }

    for (uint32_t nid = 0; nid < nid2types.size(); ++nid) {
        auto & types = nid2types[nid];
        int min_count = 0;
        string min_type; // use empty node type if no type
        for (auto & t : types) {
            uint32_t tid = type2typeId[t];
            int cur_count = typeId2Count[tid].size();
            if (min_count == 0) {
                min_count = cur_count;
                min_type = t;
                continue;
            }
            if (min_count > cur_count) {
                min_count = cur_count;
                min_type = t;
            }
        }

        nid2mintype[nid] = min_type;


    }
    auto s2 = getTime();

    printf("Time spent on finding most specific node types: %.3lf ms.\n",
           getInterval(s1, s2));

}

void GraphManager::serializeToDisk(const char *filename) {
    ofstream ofs(filename);
    if (!ofs.is_open()) {
        cout << "Invalid File path\n";
        return;
    }

    ofs.write((char *) &node_sz, sizeof(size_t));
    ofs.write((char *) &edge_sz, sizeof(size_t));

    ofs.write((char *) nodes, sizeof(uint32_t) * (node_sz + 1));
    ofs.write((char *) edges, sizeof(uint64_t) * edge_sz);
    ofs.write((char *) edgesSpreadRatio, sizeof(double) * edge_sz);


    boost::archive::binary_oarchive oa(ofs);
    oa << typeId2Name;
    oa << typeId2Count;
    oa << nid2IRI;
    oa << nid2types;
    oa << type2nid;
    oa << snidEid2pos;
    oa << iri2Nid;

    ofs.close();
}

void GraphManager::deserializeFromDisk(const char *filename) {
    ifstream ifs(filename);
    if (!ifs.is_open()) {
        printf("Error: Invalid File path For Deserialization: %s\n", filename);
        return;
    }

    ifs.read((char *) &node_sz, sizeof(size_t));
    ifs.read((char *) &edge_sz, sizeof(size_t));

    nodes = new uint32_t[node_sz + 1];
    edges = new uint64_t[edge_sz];
    edgesSpreadRatio = new double[edge_sz];

    ifs.read((char *) nodes, sizeof(uint32_t) * (node_sz + 1));
    ifs.read((char *) edges, sizeof(uint64_t) * edge_sz);
    ifs.read((char *) edgesSpreadRatio, sizeof(double) * edge_sz);


    boost::archive::binary_iarchive ia(ifs);
    ia >> typeId2Name;
    ia >> typeId2Count;
    ia >> nid2IRI;
    ia >> nid2types;
    ia >> type2nid;
    ia >> snidEid2pos;
    ia >> iri2Nid;

    ifs.close();
}


void GraphManager::enumWithRepeat(uint32_t source) {
    vector<unordered_map<string, int>> hopwiseEdge2freq((size_t) num_hop);

    vector<uint32_t> initPath;
    initPath.push_back(source);
    vector<vector<uint32_t>> pathInstances;
    vector<vector<uint32_t>> allPathInstances;
    pathInstances.push_back(initPath);
    for (int i = 0; i < num_hop; ++i) {
        vector<vector<uint32_t>> newPathInstances;
        for (int j = 0; j < pathInstances.size(); ++j) {
            auto path = pathInstances[j];
            //only consider simple paths
            unordered_set<uint32_t> visited;
            for (int k = 0; k < path.size(); ++k)
                if (k % 2 == 0)
                    visited.insert(path[k]);
            //start expansion
            uint32_t frontier = path.back();
            int adj_sz = nodes[frontier + 1] - nodes[frontier];
            for (int l = 0; l < adj_sz; ++l) {
                auto edge = edges[nodes[frontier] + l];
                uint32_t tid = extractHigh32bits(edge);
                uint32_t end_node = extractLow32bits(edge);
                if (visited.count(end_node)) continue;
                string edge_type = typeId2Name[tid];

                if (hopwiseEdge2freq[i].count(edge_type))
                    hopwiseEdge2freq[i][edge_type]++;
                else
                    hopwiseEdge2freq[i][edge_type] = 1;

                vector<uint32_t> tmp = path;
                tmp.push_back(tid);
                tmp.push_back(end_node);
                newPathInstances.push_back(tmp);
            }
        }

        pathInstances = newPathInstances;
        allPathInstances.insert(allPathInstances.end(), pathInstances.begin(), pathInstances.end());
    }

    cout << allPathInstances.size() << endl;

    //Traverse all the path instances. Let's see how pointwise mutual information works
    unordered_set<string> visited;
    for (int i = 0; i < allPathInstances.size(); ++i) {
        auto ppattern = allPathInstances[i];
        string sig = path2sig(ppattern);
//        if (visited.count(sig) || ppattern[1] != 31 || ppattern.size()!=7) continue;
        //if (visited.count(sig) || ppattern.back() != 3310) continue;
        if (visited.count(sig) || ppattern.back() != 276746) continue;
        visited.insert(sig);
        auto pathSeq = extractPathSeqOnly(ppattern);
        unordered_map<uint32_t, double> s2count;
//        singleSourcePatternCount(source, pathSeq, s2count);
        singleSourcePatternExpansionWithSink(source, pathSeq, s2count);
//        unordered_map<uint32_t, int> globalPatternCount;
//        pathPatternGlobalCount(pathSeq, globalPatternCount);


//        double p_x = 0;

//        for (auto kv_pair : s2count) {
//            p_x += kv_pair.second;
//        }
        //double ttl = 0;
//        for (auto kv_pair : globalPatternCount) {
//            ttl += kv_pair.second;
//        }
//
//        p_x /= ttl;

        //traverse all end node
        for (auto kv_pair : s2count) {

            auto sequenLabel = sig2labelSeq(sig);
            auto end_node = kv_pair.first;
            //if (end_node != 3310) continue;
            if (end_node != 276746) continue;

            unordered_map<uint32_t, double> end2count;
            auto reversePathSeq = extractPathSeqOnly(reversePathSequence(ppattern));
//            singleSourcePatternCount(end_node, reversePathSeq, end2count);
            singleSourcePatternExpansionWithSink(end_node, reversePathSeq, end2count);
            //calculate path score. Self-design
//            double p_yGivenx = kv_pair.second / p_x;

            double p_yGivenx = kv_pair.second;

//            double p_y = 0;
//            for (auto kv_pair_end : end2count) {
//                p_y += kv_pair_end.second;
//            }
            double p_xGiveny = end2count[source];
            double finalScore = (p_yGivenx + p_xGiveny) / 2;
//            if (end_node == 3310)
//            if (finalScore > 0.4)
            printf("sequenceLabel -> endNode, (pmi, npmi) = %s -> %s, (%lf, %lf, %lf)\n",
                   sequenLabel.substr(0, sequenLabel.size() - 1).c_str(), nid2IRI[end_node].c_str(), p_yGivenx,
                   p_xGiveny, finalScore);


            //calculate path pattern score.  PMI
//            double p_y = globalPatternCount[end_node];
//            p_y /= ttl;
//            double p_joint = kv_pair.second / ttl;
//
//            double pmi = log2(p_joint/(p_x * p_y));
//            double npmi = -pmi / log2(p_joint);
//
//            if (end_node == 60)
//                printf("sequenceLabel -> endNode, (pmi, npmi) = %s -> %s, (%lf, %lf)\n",
//                   sequenLabel.substr(0, sequenLabel.size() - 1).c_str(), nid2IRI[end_node].c_str(), pmi, npmi);

        }


    }

}


void GraphManager::singleSourcePatternExpansionWithSink(uint32_t source,
                                                        std::vector<uint32_t> &pathPattern,
                                                        std::unordered_map<uint32_t, double> &endNode2Prob) {
    //filter out impossible path patterns
    if (pathPattern.empty() || pathPattern[0] >= typeId2Count.size())
        return;

    endNode2Prob.clear();
//    unordered_map<uint32_t, double> frontiers;
    unordered_map<uint32_t, pair<unordered_set<uint32_t>, double>> frontiers;
    frontiers[source] = make_pair(unordered_set<uint32_t>({source}), 1);

    for (int i = 0; i < pathPattern.size(); ++i) {
        unordered_map<uint32_t, pair<unordered_set<uint32_t>, double>> tmp_frontiers;
        for (auto &kv_pair : frontiers) {
            auto f = kv_pair.first;
            double prob = kv_pair.second.second;
            int edge_start = nodes[f];
            int adj_sz = nodes[f + 1] - nodes[f];
            vector<uint32_t> possibleNeighbors;
            for (int j = 0; j < adj_sz; ++j) {
                auto edge = edges[edge_start + j];
                auto tid = extractHigh32bits(edge);
                auto end_node = extractLow32bits(edge);
                if (tid == pathPattern[i] && !kv_pair.second.first.count(end_node)) {
                    double prob2assign = prob * edgesSpreadRatio[edge_start + j];
                    if (tmp_frontiers.count(end_node)) {
                        tmp_frontiers[end_node].second += prob2assign;
                        tmp_frontiers[end_node].first.insert(f);
                    } else
                        tmp_frontiers[end_node] = make_pair(unordered_set<uint32_t>({f}), prob2assign);

                }
            }
        }

        frontiers = tmp_frontiers;
    }

    //insert last level frontiers to return
    for (auto &kv_pair : frontiers) {
        endNode2Prob[kv_pair.first] = kv_pair.second.second;
    }
//    endNode2Prob = frontiers;
}


void GraphManager::enumWithPriority(uint32_t source) {
    unordered_map<string, unordered_set<priority_entry>> patternSig2frontiers; //local frontiers from source
    priority_queue<priority_entry, std::vector<priority_entry>, pe_comparator> priorityNodePattern;


    priority_entry source_entry(source, 1, 1, "");
    priorityNodePattern.push(source_entry);
    patternSig2frontiers[""] = unordered_set<priority_entry>({source_entry});

    unordered_set<string> visited_patterns;
    unordered_map<string, pair<unordered_set<uint32_t>, int>> patternSig2occur; // GLOBAL frontiers and, global pattern count
    patternSig2occur[""] = make_pair(unordered_set<uint32_t>({source}), 1);
    int count = 0;
    while (!priorityNodePattern.empty()) {
        count++;
        auto e = priorityNodePattern.top();
        priorityNodePattern.pop();

        //can print the info, if wanted. Until output, we calculate its reverse received probability
        //lazy-calculation
        auto ppattern = sig2pathPattern(e.patternSig);
//        auto reverse_ppattern = reversePatternSequence(ppattern);
//        unordered_map<uint32_t, double> end2prob;
//        singleSourcePatternExpansionWithSink(e.nid, reverse_ppattern, end2prob);
        if (!e.patternSig.empty()) {
            //it can only be empty when expanding the source node.
//            e.score = (end2prob[source] + e.received_prob) / 2;
            printf("%s -> %s, (YFromX, XFromY, GlobalSignificance, FinalScore) = (%lf, %lf, %d, %lf)\n",
                   sig2labelSeq(e.patternSig).c_str(), nid2IRI[e.nid].c_str(),
                   e.received_prob, e.toSource_prob, e.globalPathPatternSignificance, e.score);
        }


        if (visited_patterns.count(e.patternSig) || ppattern.size() >= num_hop)
            continue;

        visited_patterns.insert(e.patternSig);

        unordered_map<uint32_t, unordered_map<uint32_t, double>> edgeLabel2NodeProb;
        unordered_map<uint32_t, unordered_map<uint32_t, unordered_set<uint32_t>>> edgeLabel2NodeParent;
        //expand the pattern
        for (auto &f_entry : patternSig2frontiers[e.patternSig]) {
            uint32_t f = f_entry.nid;
            int edge_start = nodes[f];
            int adj_sz = nodes[f + 1] - edge_start;

            //enumerate new edge pattern
            for (int m = 0; m < adj_sz; ++m) {
                auto edge = edges[edge_start + m];
                auto tid = extractHigh32bits(edge);
                if (tid % 2 != 0) continue;//only consider forward edges
                auto end_node = extractLow32bits(edge);

                if (f_entry.parents.count(end_node)) continue; //regard parents as a sink
                if (!edgeLabel2NodeProb.count(tid))
                    edgeLabel2NodeProb[tid] = unordered_map<uint32_t, double>();
                double end_node_prob = f_entry.received_prob * edgesSpreadRatio[edge_start + m];

                //insert parent
                edgeLabel2NodeParent[tid][end_node].insert(f_entry.nid);
                if (edgeLabel2NodeProb[tid].count(end_node)) {
                    edgeLabel2NodeProb[tid][end_node] += end_node_prob;
                } else {
                    edgeLabel2NodeProb[tid][end_node] = end_node_prob;
                }
            }
        }

        //push into queue and insert into 'patternSig2frontiers'
        for (auto &kv_pair : edgeLabel2NodeProb) {
            string newPathPatternSig;
            if (e.patternSig.empty()) {
                newPathPatternSig = to_string(kv_pair.first);
                patternSig2occur[e.patternSig].first = typeId2Count[kv_pair.first];
                patternSig2occur[e.patternSig].second = (int) typeId2Count[kv_pair.first].size();
            } else
                newPathPatternSig = e.patternSig + " " + to_string(kv_pair.first);


            //calculate global occurrences
            unordered_set<uint32_t> newGlobal_frontiers;
            int globalCount = globalPatternGrowCount(patternSig2occur[e.patternSig].first, newGlobal_frontiers,
                                                     kv_pair.first, patternSig2occur[e.patternSig].second);

            patternSig2occur[newPathPatternSig] = make_pair(newGlobal_frontiers, globalCount);

            for (auto &kv_pair_n2prob : kv_pair.second) {
                priority_entry tmp_pe(kv_pair_n2prob.first, kv_pair_n2prob.second, -1, newPathPatternSig,
                                      globalCount);
                //calculate score. eager-calculation
                auto tmp_ppattern = sig2pathPattern(tmp_pe.patternSig);
                auto tmp_reverse_ppattern = reversePatternSequence(tmp_ppattern);
                unordered_map<uint32_t, double> tmp_end2prob;
                singleSourcePatternExpansionWithSink(tmp_pe.nid, tmp_reverse_ppattern, tmp_end2prob);
                tmp_pe.toSource_prob = tmp_end2prob[source];
                tmp_pe.score = (tmp_pe.received_prob + tmp_pe.toSource_prob) / 2;

//                tmp_pe.score *= log2(tmp_pe.globalPathPatternSignificance);


                tmp_pe.parents.insert(edgeLabel2NodeParent[kv_pair.first][tmp_pe.nid].begin(),
                                      edgeLabel2NodeParent[kv_pair.first][tmp_pe.nid].end());//find the correct parent
                priorityNodePattern.push(tmp_pe);
                patternSig2frontiers[newPathPatternSig].insert(tmp_pe);

            }
        }
        //to save space, delete what is not used in the future
        // patternSig2frontiers.erase(e.patternSig);

        patternSig2frontiers.erase(e.patternSig);
        patternSig2occur.erase(e.patternSig);
    }

}


int GraphManager::globalPatternGrowCount(std::unordered_set<uint32_t> &frontiers,
                                         std::unordered_set<uint32_t> &frontiersAfterExpand,
                                         uint32_t eid2app,
                                         int currentCount) {
    //based on vertex overlap
    for (auto &f: frontiers) {
        int edge_start = nodes[f];
        int adj_sz = nodes[f + 1] - edge_start;
        for (int i = 0; i < adj_sz; ++i) {
            auto edge = edges[edge_start + i];
            auto tid = extractHigh32bits(edge);
            if (tid != eid2app) continue;

            frontiersAfterExpand.insert(extractLow32bits(edge));
        }
    }
    return frontiersAfterExpand.size() <= currentCount ? (int) frontiersAfterExpand.size() : currentCount;

}


void GraphManager::topkPathPatterns(uint32_t source, uint32_t target, int patternTopk, int peerNodeTopk) {
    //enumerate path patterns
    unordered_map<string, unordered_set<priority_entry>> patternSig2frontiers; //local frontiers from source

    typedef priority_queue<priority_entry, std::vector<priority_entry>, pe_comparator_pairNode> priorityQueue;
    typedef priority_queue<priority_entry, std::vector<priority_entry>, pe_comparator_pairNodeTopkHeap> topkHeap;

    priorityQueue priorityNodePattern;

    priority_entry source_entry(source, 1, 1, "");
    priorityNodePattern.push(source_entry);
    patternSig2frontiers[""] = unordered_set<priority_entry>({source_entry});

    unordered_set<string> visited_patterns;
    unordered_map<string, pair<unordered_set<uint32_t>, int>> patternSig2occur; // GLOBAL frontiers and, global pattern count
    patternSig2occur[""] = make_pair(unordered_set<uint32_t>({source}), 1);
    int count = 0;
    vector<topkHeap> topkPatterns; // to save answers
    pe_comparator_pairNodeTopkHeap comp;

    auto s_time = getTime();
    while (!priorityNodePattern.empty() && topkPatterns.size() < patternTopk) {
        count++;
        auto e = priorityNodePattern.top();
        priorityNodePattern.pop();

        //can print the info, if wanted. Until output, we calculate its reverse received probability
        //lazy-calculation
        auto ppattern = sig2pathPattern(e.patternSig);
//        auto reverse_ppattern = reversePatternSequence(ppattern);
//        unordered_map<uint32_t, double> end2prob;
//        singleSourcePatternExpansionWithSink(e.nid, reverse_ppattern, end2prob);
//        if (!e.patternSig.empty()) {
//            //it can only be empty when expanding the source node.
////            e.score = (end2prob[source] + e.received_prob) / 2;
//            printf("%s -> %s, (YFromX, XFromY, GlobalSignificance, FinalScore) = (%lf, %lf, %d, %lf)\n",
//                   sig2labelSeq(e.patternSig).c_str(), nid2IRI[e.nid].c_str(),
//                   e.received_prob, e.toSource_prob, e.globalPathPatternSignificance, e.score);
//        }
        //if (topkPatterns.size() >= patternTopk) break;
        if (visited_patterns.count(e.patternSig) || ppattern.size() >= num_hop)
            continue;

        visited_patterns.insert(e.patternSig);

        unordered_map<uint32_t, unordered_map<uint32_t, double>> edgeLabel2NodeProb;
        unordered_map<uint32_t, unordered_map<uint32_t, unordered_set<uint32_t>>> edgeLabel2NodeParent;
        //expand the pattern by adding an edge pattern
        for (auto &f_entry : patternSig2frontiers[e.patternSig]) {
            uint32_t f = f_entry.nid;
            int edge_start = nodes[f];
            int adj_sz = nodes[f + 1] - edge_start;

            //enumerate new edge pattern
            for (int m = 0; m < adj_sz; ++m) {
                auto edge = edges[edge_start + m];
                auto tid = extractHigh32bits(edge);
                //if (tid % 2 != 0) continue;//only consider forward edges
                auto end_node = extractLow32bits(edge);

                if (f_entry.parents.count(end_node)) continue; //regard parents as a sink
                if (!edgeLabel2NodeProb.count(tid))
                    edgeLabel2NodeProb[tid] = unordered_map<uint32_t, double>();
                double end_node_prob = f_entry.received_prob * edgesSpreadRatio[edge_start + m];

                //insert parent
                edgeLabel2NodeParent[tid][end_node].insert(f_entry.nid);
                if (edgeLabel2NodeProb[tid].count(end_node)) {
                    edgeLabel2NodeProb[tid][end_node] += end_node_prob;
                } else {
                    edgeLabel2NodeProb[tid][end_node] = end_node_prob;
                }
            }
        }

        //push into queue and insert into 'patternSig2frontiers'
        for (auto &kv_pair : edgeLabel2NodeProb) {
            string newPathPatternSig;
            if (e.patternSig.empty()) {
                newPathPatternSig = to_string(kv_pair.first);
                patternSig2occur[e.patternSig].first = typeId2Count[kv_pair.first];
                patternSig2occur[e.patternSig].second = (int) typeId2Count[kv_pair.first].size();
            } else
                newPathPatternSig = e.patternSig + " " + to_string(kv_pair.first);


            //calculate global occurrences
            unordered_set<uint32_t> newGlobal_frontiers;
            int globalCount = globalPatternGrowCount(patternSig2occur[e.patternSig].first, newGlobal_frontiers,
                                                     kv_pair.first, patternSig2occur[e.patternSig].second);

            patternSig2occur[newPathPatternSig] = make_pair(newGlobal_frontiers, globalCount);

            //for the same kv_pair.first that is the edge type id, the path pattern will also be the same
            //so we can calculate the peers of target nodes here.
            bool considerThisPattern = kv_pair.second.count(target);
            topkHeap tmp_queue;

            for (auto &kv_pair_n2prob : kv_pair.second) {
                priority_entry tmp_pe(kv_pair_n2prob.first, kv_pair_n2prob.second, -1, newPathPatternSig,
                                      globalCount);

                //calculate score. eager-calculation

                auto tmp_ppattern = sig2pathPattern(tmp_pe.patternSig);
                auto tmp_reverse_ppattern = reversePatternSequence(tmp_ppattern);
                unordered_map<uint32_t, double> tmp_end2prob;
                singleSourcePatternExpansionWithSink(tmp_pe.nid, tmp_reverse_ppattern, tmp_end2prob);
                tmp_pe.toSource_prob = tmp_end2prob[source];
                tmp_pe.score = (tmp_pe.received_prob + tmp_pe.toSource_prob) / 2;

//                    tmp_pe.score *= log2(tmp_pe.globalPathPatternSignificance);
                if (considerThisPattern) {
                    //insert into queue
                    if (tmp_queue.size() < peerNodeTopk)
                        tmp_queue.push(tmp_pe);
                    else {
                        if (comp(tmp_pe, tmp_queue.top())) {
                            tmp_queue.pop();
                            tmp_queue.push(tmp_pe);
                        }
                    }

                    visited_patterns.insert(tmp_pe.patternSig); // no more expansion from the already identified.
                    continue;
                }
                tmp_pe.parents.insert(edgeLabel2NodeParent[kv_pair.first][tmp_pe.nid].begin(),
                                      edgeLabel2NodeParent[kv_pair.first][tmp_pe.nid].end());//find the correct parent
                priorityNodePattern.push(tmp_pe);
                patternSig2frontiers[newPathPatternSig].insert(tmp_pe);

            }
            if (considerThisPattern) {
                topkPatterns.push_back(tmp_queue);
                if (topkPatterns.size() >= patternTopk) break; //process finish
            }
        }
        //to save space, delete what is not used in the future
        // patternSig2frontiers.erase(e.patternSig);

        patternSig2frontiers.erase(e.patternSig);
        patternSig2occur.erase(e.patternSig);

    }
    auto e_time = getTime();
    printf("Finish processing. Total time = %.3lf ms.\n", getInterval(s_time, e_time));
    //print the result
    for (int i = 0; i < topkPatterns.size(); ++i) {
        auto q = topkPatterns[i];
        printf("%d : %s \n", i, sig2labelSeq(q.top().patternSig).c_str());
        vector<priority_entry> tmp;
        while (!q.empty()) {
            auto e = q.top();
            q.pop();
            tmp.push_back(e);
        }
        reverse(tmp.begin(), tmp.end());
        for (auto &e: tmp) {
            printf("\t%s : (YFromX, XFromY, GlobalSignificance, FinalScore) = (%lf, %lf, %d, %lf)\n",
                   nid2IRI[e.nid].c_str(), e.received_prob, e.toSource_prob, e.globalPathPatternSignificance, e.score);

        }

    }
}

void
GraphManager::enumeratePathInstancesWithPPRScore(uint32_t source,
                                                 std::unordered_map<uint32_t, double> &nodePPR,
                                                 vector<vector<uint32_t>> &candidatePathIns) {
    // BFS enumeration
    vector<uint32_t> initPath;
    initPath.push_back(source);
    vector<vector<uint32_t>> pathInstances;
    vector<vector<uint32_t>> allPathInstances;
    pathInstances.push_back(initPath);

    for (int i = 0; i < num_hop; ++i) {
        vector<vector<uint32_t>> newPathInstances;
        for (int j = 0; j < pathInstances.size(); ++j) {
            auto path = pathInstances[j];
            //only consider simple paths without node repetitions
            unordered_set<uint32_t> visited;
            for (int k = 0; k < path.size(); ++k)
                if (k % 2 == 0)
                    visited.insert(path[k]);
            //start expansion
            uint32_t frontier = path.back();
            int adj_sz = nodes[frontier + 1] - nodes[frontier];
            for (int l = 0; l < adj_sz; ++l) {
                auto edge = edges[nodes[frontier] + l];
                uint32_t tid = extractHigh32bits(edge);
                uint32_t end_node = extractLow32bits(edge);
                if (visited.count(end_node)) continue;
                string edge_type = typeId2Name[tid];

                vector<uint32_t> tmp = path;//the path preserves the edge label
                tmp.push_back(tid);
                tmp.push_back(end_node);
                newPathInstances.push_back(tmp);
            }
        }

        pathInstances = newPathInstances;
        allPathInstances.insert(allPathInstances.end(), pathInstances.begin(), pathInstances.end());
    }
    vector<pair<double, int>> scoreIdx;
    for (int i = 0; i < allPathInstances.size(); ++i) {
        auto &p = allPathInstances[i];
        double avg_score = 0;
        for (int j = 0; j < p.size(); j += 2) {
            if (nodePPR.count(p[j]))
                avg_score += nodePPR[p[j]];
        }
        avg_score /= (p.size() / 2 + 1);
        scoreIdx.emplace_back(make_pair(avg_score, i));
    }

    sort(scoreIdx.begin(), scoreIdx.end());
    reverse(scoreIdx.begin(), scoreIdx.end());
    //print top-k
    int topk = 100;

    for (int m = 0; m < scoreIdx.size(); ++m) {
        //if (allPathInstances[scoreIdx[m].second].back() != 3198) continue;
        if (allPathInstances[scoreIdx[m].second].back() != 276746) continue;
        auto p = allPathInstances[scoreIdx[m].second];
        string pathLabelSeq = pathInstance2stringSeq(p);
        double s = scoreIdx[m].first;
        candidatePathIns.push_back(p);
        printf("%s = %.5lf\n", pathLabelSeq.c_str(), s);
    }

}


void GraphManager::setEdgeType2Pos() {

    ifstream ifs("snidEid2pos.txt");
    if (ifs.is_open()) {
        boost::archive::binary_iarchive ia(ifs);
        ia >> snidEid2pos;
        ifs.close();
        return;
    }
    ifs.close();
    for (uint32_t i = 0; i < node_sz; ++i) {
        int edge_start = nodes[i];
        int adj_sz = nodes[i + 1] - edge_start;
        for (int j = 0; j < adj_sz; ++j) {
            auto edge = edges[edge_start + j];
            auto eid = extractHigh32bits(edge);
            auto snidEid = combine2u32(i, eid);
            if (snidEid2pos.count(snidEid)) {
                continue;
            }
            snidEid2pos[snidEid] = edge_start + j;
        }
    }
    ofstream ofs("snidEid2pos.txt");
    boost::archive::binary_oarchive oa(ofs);
    oa << snidEid2pos;
    ofs.close();
}

void GraphManager::setIRI2NidMap() {
    ifstream ifs("iri2Nid.txt");
    if (ifs.is_open()) {
        boost::archive::binary_iarchive ia(ifs);
        ia >> iri2Nid;
        ifs.close();
        return;
    }
    ifs.close();
    string prefix = "https://www.wikidata.org/wiki/";
    for (uint32_t i = 0; i < node_sz; ++i) {
        auto iri = nid2IRI[i].substr(prefix.size() + 1);
        iri2Nid[iri] = i;
    }
    ofstream ofs("iri2Nid.txt");
    boost::archive::binary_oarchive oa(ofs);
    oa << iri2Nid;
    ofs.close();


}

void GraphManager::sortSameTypeNodesByDegree() {
    for (int i = 0; i < node_sz; ++i) {
        int edge_start = nodes[i];
        int adj_sz = nodes[i + 1] - nodes[i];
        int iter_begin = edge_start;
        int iter_end = edge_start;
        vector<pair<int, uint64_t>> sameTypedEdge;
        for (int j = 0; j < adj_sz; ++j) {
            auto edg = edges[edge_start + j];
            auto eid = extractHigh32bits(edg);
            auto end_node = extractLow32bits(edg);
            int end_node_degree = nodes[end_node + 1] - nodes[end_node];
            if (j == 0 || eid == extractHigh32bits(edges[edge_start + j - 1])) {
                sameTypedEdge.emplace_back(make_pair(end_node_degree, edg));
                continue;
            }
            sort(sameTypedEdge.begin(), sameTypedEdge.end(), greater<pair<int, uint64_t>>());
            //copy back
            for (int k = 0; k < sameTypedEdge.size(); ++k) {
                int idx_1 = (int)sameTypedEdge.size() - 1 - k;
                int idx_2 = j - 1 -k;
                edges[edge_start + idx_2] = sameTypedEdge[idx_1].second;
            }
            sameTypedEdge.clear();
            sameTypedEdge.emplace_back(make_pair(end_node_degree, edg));
        }

    }


}