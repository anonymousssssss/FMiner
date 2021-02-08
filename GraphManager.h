//
// Created by yangyueji on 7/15/19.
//

#ifndef PATHPATTERNMINING_GRAPHMANAGER_H
#define PATHPATTERNMINING_GRAPHMANAGER_H


#include <cstdio>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <sys/time.h>


inline
unsigned long long getTime() {
    struct timeval tv;

    gettimeofday(&tv, NULL);

    unsigned long long ret = tv.tv_usec;

    /* Adds the seconds after converting them to microseconds (10^-6) */
    ret += (tv.tv_sec * 1000 * 1000);

    return ret;
};



inline
double getInterval(unsigned long long start, unsigned long long stop) {
    // return millisecond
    return (double) (stop - start) / 1000.0;
};


__inline__
uint32_t extractHigh32bits(uint64_t num) {
    return (uint32_t) ((num >> 32) & 0xFFFFFFFF);
}

__inline__
uint32_t extractLow32bits(uint64_t num) {
    return (uint32_t) (num & 0xFFFFFFFF);
}

__inline__
uint64_t combine2u32(uint32_t high, uint32_t low) {
    uint64_t res = 0;
    res = res | high;
    res = res << 32;
    res = res | low;
    return res;
}

//template <typename T1, typename T2, typename T3>
struct priority_entry {
    uint32_t nid;
    double received_prob; //probability received from source
    double toSource_prob;
    double score;// probability assigned. Bidirectional
    std::string patternSig;
    int globalPathPatternSignificance;
    std::unordered_set<uint32_t> parents; //to avoid going immediately back


    priority_entry(uint32_t n, double rp, double s, std::string p, int rpc = 1, double toSourceP = 0) {
        nid = n;
        received_prob = rp;
        toSource_prob = toSourceP;
        score = s; // -1 means not calculated yet
        patternSig = p;
        globalPathPatternSignificance = rpc;
    }

    bool operator==(const priority_entry &pe2) const {
        return nid == pe2.nid && patternSig == pe2.patternSig;
    }
};


namespace std {
    template<>
    struct hash<priority_entry> {
        typedef priority_entry argument_type;
        typedef std::size_t result_type;

        result_type operator()(argument_type const &pe) const noexcept {
            result_type const h(std::hash<uint32_t>{}(pe.nid));
            return h;
        }
    };
}

class pe_comparator {
public:
    bool operator()(const priority_entry &p1, const priority_entry &p2) {
        /// based on total score
//        return p1.score < p2.score;
        ///based on received prob
//        return p1.received_prob < p2.received_prob;
        ///linear scoring on score and significance.
        if (p1.score < p2.score) return true;
        if (p1.score > p2.score) return false;
        return p1.globalPathPatternSignificance < p2.globalPathPatternSignificance;

    }
};

class pe_comparator_pairNode {
public:
    bool operator()(const priority_entry &p1, const priority_entry &p2) {
        /// based on total score
//        return p1.score < p2.score;
        ///based on received prob
//        return p1.received_prob < p2.received_prob;
        ///linear scoring on score and significance.
        if (p1.score < p2.score) return true;
        if (p1.score > p2.score) return false;
        return p1.globalPathPatternSignificance < p2.globalPathPatternSignificance;

    }
};

class pe_comparator_pairNodeTopkHeap {
public:
    bool operator()(const priority_entry &p1, const priority_entry &p2) {
        /// based on total score
//        return p1.score > p2.score;
        ///based on received prob
//        return p1.received_prob > p2.received_prob;
        ///linear scoring on score and significance.
        if (p1.score > p2.score) return true;
        if (p1.score < p2.score) return false;
        return p1.globalPathPatternSignificance > p2.globalPathPatternSignificance;

    }
};

class GraphManager {
public:
    /**Graph Adjacency Lists.
     * Type specific store.*/
    size_t node_sz = 0; // 1 larger than the actual node number
    size_t edge_sz = 0;

    uint32_t *nodes = nullptr;
    uint64_t *edges = nullptr; //high 32-bit for edge type, low 32-bit for edge id
    double *edgesSpreadRatio = nullptr; //w.r.t. a specific type of edge


    /***Node type maps*/
//    std::vector<std::vector<std::string>> nid2types;

    std::vector<std::unordered_set<std::string>> nid2types;
    std::unordered_map<std::string, std::unordered_set<uint32_t>> type2nid;
    std::vector<std::string> nid2IRI;
    std::unordered_map<uint32_t, std::string> nid2mintype;
    /**Edge Type to start node.*/
//    size_t type_sz = 0;
    std::vector<std::string> typeId2Name;
    std::vector<std::unordered_set<uint32_t>> typeId2Count; // edge type to source node

    int num_hop = 3;

    GraphManager() = default;

    ~GraphManager() {
        delete[] nodes;
        delete[] edges;
        delete[] edgesSpreadRatio;
    }

    //load the graph from edge file
    void readEdges(const char *filename);

    //load node types
    void readNodeTypes(const char * filename);
    void resetNodeTypes();

    // serialization implementation
    void serializeToDisk(const char *filename);

    void deserializeFromDisk(const char *filename);

    void enumWithRepeat(uint32_t source);

    //enumerate all path patterns and nodes neighboring the source
    void enumWithPriority(uint32_t source);

    //directly pass in the path sequence without considering nodes.
    void singleSourcePatternExpansionWithSink(uint32_t source, std::vector<uint32_t> &pathPattern,
                                              std::unordered_map<uint32_t, double> &endNode2Prob);

    //count global significant for a given pattern and an edge to expand.
    int globalPatternGrowCount(std::unordered_set<uint32_t> &frontiers,
                               std::unordered_set<uint32_t> &frontiersAfterExpand,
                               uint32_t eid2app,
                               int currentCount);

    /**next to implement two functions:
     * 1. Given a pair of nodes, find the top-k path patterns.
     * 2. Given the top-k path patterns find the top-k similar nodes within each path patterns.
     * The second step can be embedded in the first step.*/
    void topkPathPatterns(uint32_t source, uint32_t target,
                          int patternTopk = 10,
                          int peerNodeTopk = 10);


    //enumerate path in a graph
    void enumeratePathInstancesWithPPRScore(uint32_t source,
                                            std::unordered_map<uint32_t, double> &nodePPR,
                                            std::vector<std::vector<uint32_t>> &candidatePathIns);

    /**Given an edge type, find the first edge position*/
    int firstEdgePos(uint32_t start_node, uint32_t eid) const;

    /** Set edge type to edge pos */
    std::unordered_map<uint64_t, int> snidEid2pos;
    void setEdgeType2Pos();

    std::unordered_map<std::string, uint32_t> iri2Nid;
    void setIRI2NidMap();

    void sortSameTypeNodesByDegree();

    inline
    std::string path2sig(std::vector<uint32_t> &path) {
        std::string line;
        for (int i = 0; i < path.size(); ++i) {
            if (i % 2 == 1) {
//                line += typeId2Name[path[i]] + ",";
                line += std::to_string(path[i]) + " ";
            }
        }
        return line;
    }

    inline
    std::string pathPattern2sig(std::vector<uint32_t> &pathPattern) {
        std::string line;
        for (int i = 0; i < pathPattern.size(); ++i) {
//                line += typeId2Name[path[i]] + ",";
            line += std::to_string(pathPattern[i]) + " ";

        }
        return line;
    }


    inline
    std::vector<uint32_t> sig2pathPattern(const std::string &sig) {
        std::vector<uint32_t> resVec;
        if (sig.empty())
            return resVec;
        uint32_t token;
        std::istringstream iss(sig);
        while (iss >> token) {
            resVec.push_back(token);
        }
        return resVec;
    }

    inline
    std::string sig2labelSeq(const std::string &sig) {
        std::istringstream iss(sig);
        std::string line;
        int eid;
        while (iss >> eid) {
            line += typeId2Name[eid] + ",";
        }
        if (!line.empty() && line.back() == ',')
            line = line.substr(0, line.length() - 1);
        return line;
    }

    inline
    uint32_t reverseEdgeTypeId(uint32_t eid) const {
        if (eid % 2 == 0)
            return eid + 1;
        else
            return eid - 1;
    }

    inline
    bool equal_pathPattern(std::vector<uint32_t> &v1, std::vector<uint32_t> &v2) {
        if (v1.size() != v2.size())
            return false;
        for (int i = 0; i < v1.size(); ++i) {
            if (i % 2 == 0) continue;
            if (v1[i] != v2[i]) return false;
        }
        return true;
    }

    inline
    std::vector<uint32_t> reversePathSequence(std::vector<uint32_t> &pathPattern) {
        std::vector<uint32_t> outPut;
        for (int i = 0; i < pathPattern.size(); ++i) {
            if (i % 2 == 1)
                outPut.push_back(reverseEdgeTypeId(pathPattern[i]));
            else
                outPut.push_back(pathPattern[i]);

        }

        std::reverse(outPut.begin(), outPut.end());
        return outPut;
    }

    inline
    std::vector<uint32_t> reversePatternSequence(std::vector<uint32_t> &pathPattern) {
        std::vector<uint32_t> outPut;
        for (int i = 0; i < pathPattern.size(); ++i)
            outPut.push_back(reverseEdgeTypeId(pathPattern[pathPattern.size()- i - 1]));
        //std::reverse(outPut.begin(), outPut.end());
        return outPut;
    }


    inline
    std::vector<uint32_t> extractPathSeqOnly(const std::vector<uint32_t> &pathPattern) {
        std::vector<uint32_t> outPut;
        for (int i = 0; i < pathPattern.size(); ++i) {
            if (i % 2 == 1)
                outPut.push_back(pathPattern[i]);
        }
        return outPut;
    }

    inline
    std::string pathInstance2stringSeq(const std::vector<uint32_t> &pathInstance) {
        std::string res;
        for (int i = 0; i < pathInstance.size(); ++i) {
            if (i % 2 == 0) //it is a node
                res += nid2IRI[pathInstance[i]] + " ,";
            else
                res += typeId2Name[pathInstance[i]] + " ,";
        }
        res = res.substr(0, res.length() - 1);
        return res;
    }

    inline void getOneLabelRandom(uint32_t nid, std::string &_label) {
        if (nid >= node_sz) return;
        for (auto & _l : nid2types[nid]) {
            _label = _l;
            return;
        }
        _label = "ANY";
    }
};


#endif //PATHPATTERNMINING_GRAPHMANAGER_H
