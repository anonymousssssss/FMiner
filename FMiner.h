//
// Created by Yueji YANG on 25/12/20.
//

#ifndef PATHPATTERNMINING_FMINER_H
#define PATHPATTERNMINING_FMINER_H

#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <cmath>
#include <utility>
#include <cassert>
#include "GraphManager.h"

// 3 minute
#define TIME_LIMIT_MS 180000

// <path score, path>, where path = (node, edge, node, edge, node....)
// The path always ends with a node id.
//typedef std::pair<double, std::vector<uint32_t>> BFS_entry;
struct BFS_entry {
    // The score is path score during path enumeration.
    // It becomes pattern score during pattern evaluation.
    double score;
    std::vector<uint32_t> path;

    BFS_entry() = default;

    BFS_entry(double score, std::vector<uint32_t> path) {
        this->score = score;
        this->path = std::move(path);
    }

    BFS_entry(double score, uint32_t init_node) {
        this->score = score;
        this->path.push_back(init_node);
    }

    bool operator==(const BFS_entry &other) const {
        bool equal = true;
        if (this->path.size() != other.path.size()) return false;
        for (int i = 0; i < this->path.size(); ++i) {
            if (this->path[i] != other.path[i])
                return false;
        }
        return equal;
    }

};

struct BFSOrder {
    bool operator()(const BFS_entry &lhs, const BFS_entry &rhs) {
        return lhs.score < rhs.score;
    }
};


// Mains a hash map used as memo.
struct pattern_memo {

    typedef std::vector<std::unordered_set<uint32_t>> matching_nodes_vec;

    std::unordered_map<std::string, std::pair<matching_nodes_vec, bool>> pattern_sig_matching_nodes;


    pattern_memo() = default;

    static
    void pattern_sig(std::vector<uint32_t> &path,
                     std::vector<std::string> &node_pattern,
                     uint32_t frontier, int path_idx,
                     std::string &pattern_sig) {
        for (int i = 0; i < path.size() && i <= path_idx; ++i) {
            if (i + 1 == path.size()) {
                pattern_sig.append(std::to_string(frontier));
            } else if (i % 2 == 0) {
                // node
                pattern_sig.append(node_pattern[i / 2]);
                pattern_sig.append(",");
            } else {
                // edge, simplify and makes results more diverse
                auto eid = path[i];
                if (path[i] % 2 > 0) {
                    eid = path[i] - 1; // uniform egdes. Fix to even numbers. Consider as undirected graphs.
                }
                pattern_sig.append(std::to_string(eid));
                pattern_sig.append(",");
            }
        }
    }

    static
    void node_pattern_sig(std::vector<uint32_t> &path,
                          std::vector<std::string> &node_pattern,
                          uint32_t frontier, int path_idx,
                          std::string &pattern_sig) {
        for (int i = 0; i < path.size() && i <= path_idx; ++i) {
            if (i + 1 == path.size()) {
                pattern_sig.append(std::to_string(frontier));
            } else if (i % 2 == 0) {
                // node
                pattern_sig.append(node_pattern[i / 2]);
                pattern_sig.append(",");
            } else {
                // edge, simplify and makes results more diverse
                continue;
            }
        }
    }


    bool has_key(std::string &_key) {
        return pattern_sig_matching_nodes.count(_key);
    }

    matching_nodes_vec &get_by_key(std::string &_key) {
        assert(has_key(_key));
        return pattern_sig_matching_nodes[_key].first;
    }

    void insert(std::string &_key, matching_nodes_vec &inp, bool exceed = false) {
        assert(!has_key(_key));
        pattern_sig_matching_nodes[_key] = std::make_pair<>(std::move(inp), exceed);
    }

    void insert(std::string &_key) {
        assert(!has_key(_key));
        // insert empty
        pattern_sig_matching_nodes[_key] = std::make_pair<>(matching_nodes_vec(), false);
    }

    bool is_exceed(std::string &_key) {
        assert(has_key(_key));
        return pattern_sig_matching_nodes[_key].second;
    }


    typedef std::vector<std::pair<std::unordered_set<uint32_t>, double>> matching_nodes_pair_vec;

    void insert(std::string &_key, matching_nodes_pair_vec &inp, bool exceed = false) {
        assert(!has_key(_key));
        pattern_sig_matching_nodes[_key].first.resize(inp.size());
        for (int i = 0; i < inp.size(); ++i) {
            pattern_sig_matching_nodes[_key].first[i] = std::move(inp[i].first);
        }
        pattern_sig_matching_nodes[_key].second = exceed;
    }

};



// Make BFS_entry hashable
namespace std {
    template<>
    struct hash<BFS_entry> {
        typedef BFS_entry argument_type;
        typedef std::size_t result_type;

        result_type operator()(argument_type const &in) const {
            size_t size = in.path.size();
            size_t seed = 0;
            for (size_t i = 0; i < size; i++)
                //Combine the hash of the current vector with the hashes of the previous ones
                //boost::hash_combine(seed, in.path[i]); Source code copied below. Otherwise need to link boost library
                seed ^= hash<uint32_t>()(in.path[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            return seed;
        }
    };
}


typedef std::priority_queue<BFS_entry, std::vector<BFS_entry>, BFSOrder> BFS_queue;

// use pointer to avoid copy and reuse the allocated memory
typedef std::unordered_map<uint32_t, std::unordered_set<BFS_entry>> node_path;

// used for storing top-k results with COFs
struct cof_entry {
    double striking_score;
    int freq_count;
    uint32_t attribute_id; // edge
    uint32_t attribute_val; // end_node

    // attribute_val with higher (or equal) frequencies.
    std::vector<std::pair<uint32_t, int>> explanations;

    double highest_freq = 0;
    std::string highest_IRI = "";

    cof_entry() = default;

    cof_entry(double score_, int freq_, uint32_t id_, uint32_t val_, std::vector<std::pair<uint32_t, int>> &exp_) :
            striking_score(score_),
            freq_count(freq_),
            attribute_id(id_),
            attribute_val(val_),
            explanations(std::move(exp_)) {};

};

struct COFOrder {
    bool operator()(const cof_entry &lhs, const cof_entry &rhs) {
        return lhs.striking_score > rhs.striking_score;
    }
};


typedef std::priority_queue<cof_entry, std::vector<cof_entry>, COFOrder> topk_cof_queue;


struct topk_entry {
    BFS_entry bfs_entry;
    std::vector<std::string> pattern;
    std::unordered_set<uint32_t> peer_entities;

    // Ranked top-k striking_score
    // <striking score, <attribute_edge, value_node>>
    std::vector<cof_entry> top_k_cofs;

    topk_entry() = default;

    // Construct without BFS_entry
    topk_entry(double score, std::vector<uint32_t> &path, std::vector<std::string> &pattern_,
               std::unordered_set<uint32_t> peer_entities_) :
            bfs_entry(std::move(BFS_entry(score, path))),
            pattern(std::move(pattern_)),
            peer_entities(std::move(peer_entities_)) {};


    void pprint(GraphManager &gm,
                std::string &path_str,
                std::string &pattern_str);

    void pprint_cof(GraphManager &gm, std::string &cof_str,
                    int target,
                    int print_num = 3,
                    int print_exp = 2);
};

struct TopkOrder {
    bool operator()(const topk_entry &lhs, const topk_entry &rhs) {
        return lhs.bfs_entry.score > rhs.bfs_entry.score;
    }
};

typedef std::priority_queue<topk_entry, std::vector<topk_entry>, TopkOrder> topk_pattern_queue;


// Used for recording visited patterns. Patterns may contain instance nodes.
// node => patterns
typedef std::unordered_map<uint32_t, std::unordered_set<std::string>> node_patterns;
// pattern (from target towards context) => matching node sets
typedef std::unordered_map<std::string,
        std::vector<std::unordered_set<uint32_t>>> pattern_matching_nodes;


class FMiner {
/***
 * Optimizations
 * 1. Bidirectional BFS.
 * 2. Delayed expansion of high-degree nodes.
 * 3. Top-k pruning.
 * 4. More accurate node type mapping.
 * 5. More accurate COF ranking.
 *    Maybe only need one for one pattern to keep diversity.
 * **/

/***
 * Notes:
 * 1. no cycle should be considered.
 * 2. assume two nodes only have one edge.
 * **/

public:
    int num_peer_entities = 2;
    int path_length_limit = 4;
    int pattern_topk_num = 20;
    int cof_topk_num = 10;
    int degree_threshold = 10000;

    // Delay high degree nodes from expansion
    bool delay_highdegree_expansion = true;

    double decay_factor = std::exp(1);

    // Set 0.1 and 2 for development purpose on sample graphs.
    // 0.95 and 500 should be used on real graph.
    // These parameter needs to be tuned.
    double attribute_occurrence_ratio_limit = 0.80;
    uint32_t degree_importance_limit = 100;

    ///For ease of access
    std::unordered_map<uint32_t, double> node_score;
    std::unordered_set<uint32_t> context_nodes;

    /// Profiling (in millisecond)
    double BFS_expansion_time = 0;

    double pattern_evaluation_time = 0;
    double pattern_eval_time2 = 0;
    int pattern_evaluation_num = 0;

    double pattern_eval_path_score_time = 0;
    double pattern_eval_node_type_time = 0;
    double pattern_enumeration_time = 0;
    double matching_node_set_call_time = 0;

    double COF_extraction_time = 0;
    int inserted_pattern_num = 0;
    int uninserted_pattern_num = 0;

    // Optimization controller
    bool use_memo = true;
    int memo_use_num = 0;

    // top-k pruned path
    int visited_path_num = 0;
    int top_k_pruned_path_num = 0;

    int not_enough_peers = 0;

    unsigned long long time_out_start;

    bool verbose = false;

    bool use_size_bound = false;
    bool size_bound_opt_on = true;
private:
    std::unordered_map<uint32_t, std::string> most_rel_node_type;
    std::unordered_map<std::string, double> type_rel_score;
    topk_pattern_queue topk_patterns;
    std::vector<topk_entry> final_topk_pattern_entries;

    // -1 means un-initialized.
    int fwd_path_length_limit = -1;
    int bwd_path_length_limit = -1;

    std::unordered_set<std::string> processed_paths;
    std::unordered_set<std::string> processed_patterns;

    pattern_memo dp_pattern_memo;

    // necessary, e.g., (14936985, 679)
    std::unordered_map<uint32_t, std::unordered_map<std::string, double>> edge_node_score_sum;

    //visited set number
    std::unordered_set<size_t> visited_peer_set;

public:
    void FMiner_biBFS(GraphManager &gm, uint32_t target);

private:

    inline
    bool is_timeout() {
        return getInterval(this->time_out_start, getTime()) >= TIME_LIMIT_MS;
    }


    inline
    void _init() {
        most_rel_node_type.clear();
        type_rel_score.clear();
        while (!topk_patterns.empty())
            topk_patterns.pop();
        final_topk_pattern_entries.clear();
        dp_pattern_memo.pattern_sig_matching_nodes.clear();

        fwd_path_length_limit = path_length_limit / 2 + 1;
        bwd_path_length_limit = (path_length_limit + 1) / 2;
        time_out_start = getTime();
    }

    double topk_insertion_time = 0;
    double size_limit_time = 0;

    double last_level_matching_time = 0;
    inline
    void on_exit() {
        if (verbose) {
            printf(">>>>> Inner PE time = %.3lf ms. \n"
                   ">>>>> Topk insertion time = %.3lf.\n"
                   ">>>>> Size limit time = %.3lf ms. \n"
                   ">>>>> Memo used times = %d",
                   pattern_eval_time2, topk_insertion_time,
                   size_limit_time, memo_use_num);
            printf(">>>>> Not enough peers=%d\n", not_enough_peers);

            /// For csv usage:
        }
        pattern_evaluation_time -= pattern_eval_node_type_time;
        if (is_timeout()) {
            pattern_evaluation_time = std::min((double) TIME_LIMIT_MS, pattern_evaluation_time);
            BFS_expansion_time = std::min((double) TIME_LIMIT_MS, BFS_expansion_time);
            COF_extraction_time = std::min((double) TIME_LIMIT_MS, COF_extraction_time);
        }

    }

    friend void topk_entry::pprint(GraphManager &gm,
                                   std::string &path_str,
                                   std::string &pattern_str);

    inline
    bool is_path_visited(std::vector<uint32_t> &path, bool as_visited = true) {
        std::string sig;
        for (int i = 0; i < path.size(); ++i) {
            sig.append(std::to_string(path[i]));
            if (i + 1 != path.size()) {
                sig.append("-");
            }
        }
        if (processed_paths.count(sig)) return true;
        if (as_visited) processed_paths.insert(sig);
        return false;
    }


    // Return true when insertion succeeds. Return false otherwise.
    inline
    bool insert_topk_pattern(double pattern_relevant_score,
                             std::vector<uint32_t> &path,
                             std::vector<std::string> &per_pattern,
                             std::unordered_set<uint32_t> &peers) {
        if (topk_patterns.size() < pattern_topk_num) {
            topk_entry tmp_topk_entry(pattern_relevant_score, path,
                                      per_pattern, peers);
            topk_patterns.emplace(tmp_topk_entry);
            inserted_pattern_num++;
            if (verbose)
                printf("######### Inserted due to not enough topk\n");
            return true;
        } else {
            if (topk_patterns.top().bfs_entry.score < pattern_relevant_score) {
                topk_entry tmp_topk_entry(pattern_relevant_score, path,
                                          per_pattern, peers);
                if (verbose)
                    printf("######### Inserted due to too small kth value %.8lf,%.8lf\n",
                           topk_patterns.top().bfs_entry.score, pattern_relevant_score);
                topk_patterns.pop();
                topk_patterns.emplace(tmp_topk_entry);
                inserted_pattern_num++;
                return true;
            }
        }
        if (verbose)
            printf("######## Insert Failed. Worse the top-k.\n");
        uninserted_pattern_num++;
        return false;
    }

    inline
    bool above_k_th_pattern_score(double score) {
        return topk_patterns.size() < pattern_topk_num || topk_patterns.top().bfs_entry.score < score;
    }

    inline
    double get_node_score(uint32_t node) {
        return node_score.count(node) ? node_score[node] : 0;
    }

    //internal usage
    void BFS_expansion(GraphManager &gm, const BFS_entry &frontier,
                       BFS_queue &bfs_queue, node_path &processed_path,
                       uint32_t target, uint32_t context,
                       BFS_queue &delayed_queue,
                       bool is_fwd = true);

    void BiBFS_expansion(GraphManager &gm,
                         BFS_queue &fwd_queue,
                         BFS_queue &bwd_queue,
                         node_path &fwd_processed_path,
                         node_path &bwd_processed_path,
                         uint32_t target, uint32_t context,
                         BFS_queue &fwd_delayed_queue,
                         BFS_queue &bwd_delayed_queue);

    static inline
    bool check_for_complete_path(uint32_t end_node, const node_path &processed_path) {
        return processed_path.count(end_node);
    }

    static inline
    bool has_cycle(std::vector<uint32_t> &path) {
        std::unordered_set<uint32_t> tmp_set;
        for (int i = 0; i < path.size(); i += 2) {
            if (tmp_set.count(path[i])) return true;
            tmp_set.insert(path[i]);
        }
        return false;
    }

    static inline
    bool node_has_type(GraphManager &gm,
                       uint32_t node,
                       std::string &type) {
        return gm.nid2types[node].count(type);
    }

    static inline
    int node_type_id(std::string &type) {
        char *p;
        int converted = (int) strtol(type.c_str(), &p, 10);
        if (*p) return -1; // it is a string type
        else return converted; // it is a number
    }

    void pattern_enumeration(std::vector<uint32_t> &path,
                             std::vector<std::string> &path_pattern,
                             size_t idx,
                             std::vector<std::vector<std::string>> &enumerated_path_pattern);

    // only need to consider this when hitting the context node.
    bool get_node_type(GraphManager &gm,
                       uint32_t target,
                       uint32_t node,
                       std::string &most_rel_type);


    void pattern_evaluation(GraphManager &gm,
                            uint32_t target,
                            const BFS_entry &entry,
                            node_path &processed_path,
                            bool is_fwd = true);


    // Find the matching node sets.
    // We allow internal nodes to be types or instance node.
    // <end_node, node_sequence> => determined_matching_node_sets
    bool matching_node_set(GraphManager &gm,
                           std::vector<uint32_t> &path,
                           std::vector<std::string> &pattern,
                           uint32_t frontier,
                           std::vector<std::pair<std::unordered_set<uint32_t>, double>> &matching_nodes,
                           int idx,
                           bool &size_exceed);



    // Extract the top-k COFs
    void topk_cof_extraction(GraphManager &gm, uint32_t target);

    void cof_scoring(GraphManager &gm, uint32_t target, topk_entry &entry);

    // Select the edges worth exploring for OFs.
    // Customize this function to achieve different strategies
    bool validate_edge(GraphManager &gm, std::unordered_set<uint32_t> &peers, uint32_t edge_tid, uint32_t end_node);


    inline
    double get_score_sum_out_edge(GraphManager &gm,
                                  uint32_t edge_type_id,
                                  std::string &node_type) {
        double score_sum = 0;
        auto &start_nodes = gm.typeId2Count[edge_type_id];

        if (edge_node_score_sum.count(edge_type_id)
            && edge_node_score_sum[edge_type_id].count(node_type)) {
            score_sum = edge_node_score_sum[edge_type_id][node_type];
            if (verbose)
                printf("stored score %ld, ", start_nodes.size());
        } else {
            // first-time calculate
            if (verbose)
                printf("%ld, ", start_nodes.size());
            for (auto &n : start_nodes) {
                // avoid score from context node.
                if (!context_nodes.count(n) && gm.nid2types[n].count(node_type)) {
                    score_sum += (node_score.count(n) ? node_score[n] : 0);
                }
            }
            edge_node_score_sum[edge_type_id][node_type] = score_sum;
        }
        return score_sum;
    }

    inline
    bool exceed_size_limit(GraphManager &gm,
                           int node_idx,
                           std::string &node_type,
                           size_t collected_size,
                           uint32_t right_edge_type,
                           int &limit_size) {
        if (collected_size > degree_threshold) {
            use_size_bound = true;
        } else {
            return false;
        }
        if (topk_patterns.size() < pattern_topk_num) {
            return false;
        }
        double k_the_score = topk_patterns.top().bfs_entry.score;
        uint32_t reversed_right_edge_type = gm.reverseEdgeTypeId(right_edge_type);
        auto sizelimit_time1 = getTime();
        double score_sum = get_score_sum_out_edge(gm, reversed_right_edge_type, node_type);
        auto sizelimit_time2 = getTime();
        size_limit_time += getInterval(sizelimit_time1, sizelimit_time2);
        double _decayed = pow(decay_factor, -node_idx);
        double sz_limit = _decayed * score_sum / k_the_score;
        limit_size = (int) sz_limit;
        return collected_size > (size_t) sz_limit;
    }
};


#endif //PATHPATTERNMINING_FMINER_H
