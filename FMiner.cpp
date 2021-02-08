//
// Created by Yueji YANG on 25/12/20.
//

#include "FMiner.h"
#include <algorithm>
#include <iostream>

using namespace std;

void
topk_entry::pprint(GraphManager &gm,
                   std::string &path_str,
                   std::string &pattern_str) {
    path_str.clear();
    pattern_str.clear();
    bool flip = true;
    for (int i = 0; i < bfs_entry.path.size(); ++i) {
        auto _id = bfs_entry.path[i];
        if (i % 2 == 0) {
            // node
            int instance_node = FMiner::node_type_id(pattern[i / 2]);

            if (flip) {
                path_str.append(gm.nid2IRI[_id]);
                if (instance_node == -1) {
                    pattern_str.append(pattern[i / 2]);
                } else {
                    pattern_str.append(gm.nid2IRI[_id]);
                }
                flip = !flip;
            } else {
                path_str.append(gm.nid2IRI[_id]);
                if (instance_node == -1) {
                    pattern_str.append(pattern[i / 2]);
                } else {
                    pattern_str.append(gm.nid2IRI[_id]);
                }
                flip = !flip;
            }
        } else {
            // edge
            path_str.append(" --(" + gm.typeId2Name[_id] + ")-->> ");
            pattern_str.append(" --(" + gm.typeId2Name[_id] + ")-->> ");
        }
    }

    pattern_str.append(".\nPeer Entities Num=" + to_string(this->peer_entities.size()));
    pattern_str.append(", Relevance Score=" + to_string(this->bfs_entry.score));
}

void
topk_entry::pprint_cof(GraphManager &gm,
                       std::string &cof_str,
                       int target,
                       int print_num,
                       int print_exp) {
    cof_str.clear();
    auto target_name = gm.nid2IRI[target];
    for (int i = 0; i < top_k_cofs.size(); ++i) {
        auto &cof = top_k_cofs[i];
        double striking_score = cof.striking_score;
        string attr_name = gm.typeId2Name[cof.attribute_id];
        string attr_val = gm.nid2IRI[cof.attribute_val];
        cof_str.append("\t\tTop-" + to_string(i + 1) + " COF: ");
        cof_str.append(target_name + " --(" + attr_name + ")-->> " + attr_val + ", striking_score = " +
                       to_string(striking_score) + "\n");
        cof_str.append("\t\t\t Explanations: ");
        int highest_freq = (int) cof.highest_freq;
        string ax = "currentFreq = " + to_string(cof.freq_count) + ", highestFreq = " + cof.highest_IRI + " : " +to_string(highest_freq);
        cof_str.append(ax);
        //        for (int j = 0; j < cof.explanations.size() && j < print_exp; ++j) {
//            auto &val_freq = cof.explanations[j];
//            cof_str.append("(" + gm.nid2IRI[val_freq.first] + ", " + to_string(val_freq.second) + "), ");
//        }
        cof_str.append("\n");
    }
}


void
FMiner::FMiner_biBFS(GraphManager &gm, uint32_t target) {
//    num_peer_entities = 2; // For debug usage.
//    degree_threshold = 400;// For debug usage.
//    verbose = true;
//    delay_highdegree_expansion= false;
//    use_memo = true;
//    size_bound_opt_on = false;

    _init();
    if (!node_score.count(target)) {
        printf("Input entities are irrelevant.");
        return;
    }

    uint32_t context;
    // For now, we only handle one context entity.
    assert(context_nodes.size() == 1);
    for (auto n : context_nodes) {
        context = n;
    }

    // forward queue
    BFS_queue fwd_queue;
    node_path fwd_processed_path; // Forward processed node path.

    // forward delayed queue
    BFS_queue fwd_delayed_queue;

    // backward queue
    BFS_queue bwd_queue;
    node_path bwd_processed_path; // Backward processed path.

    // forward delayed queue
    BFS_queue bwd_delayed_queue;


    fwd_queue.push(BFS_entry(node_score[target], target));
    bwd_queue.push(BFS_entry(node_score[context], context));

    BiBFS_expansion(gm,
                    fwd_queue,
                    bwd_queue,
                    fwd_processed_path,
                    bwd_processed_path,
                    target, context,
                    fwd_delayed_queue,
                    bwd_delayed_queue);

    if (this->delay_highdegree_expansion && !is_timeout()) {
        this->delay_highdegree_expansion = false;
        BiBFS_expansion(gm,
                        fwd_delayed_queue,
                        bwd_delayed_queue,
                        fwd_processed_path,
                        bwd_processed_path,
                        target, context,
                        fwd_delayed_queue,
                        bwd_delayed_queue);
    }

    // path and pattern enumeration is done. Now, use top-k patterns to generate top-k COFs.
    auto topk_cof1 = getTime();
    topk_cof_extraction(gm, target);
    auto topk_cof2 = getTime();
    COF_extraction_time += getInterval(topk_cof1, topk_cof2);
//
//    // print the result. Only print 3 patterns.
    for (int i = 0; i < final_topk_pattern_entries.size(); ++i) {
        auto &e = final_topk_pattern_entries[i];
        string path_str, pattern_str, fact_str;
        e.pprint(gm, path_str, pattern_str);
        e.pprint_cof(gm, fact_str, target);
        printf("\nTop-%d relevant pattern: \n\tpath: %s, \n\tpattern: %s\n", i, path_str.c_str(), pattern_str.c_str());
        printf("%s", fact_str.c_str());
    }
//    printf("\n\tForward visited path = %ld, bacward visited path = %ld, Memo used numbers = %d.", fwd_processed_path.size(), bwd_processed_path.size(), memo_use_num);
    on_exit();
}

void
FMiner::BiBFS_expansion(GraphManager &gm,
                        BFS_queue &fwd_queue,
                        BFS_queue &bwd_queue,
                        node_path &fwd_processed_path,
                        node_path &bwd_processed_path,
                        uint32_t target, uint32_t context,
                        BFS_queue &fwd_delayed_queue,
                        BFS_queue &bwd_delayed_queue) {
    while (!fwd_queue.empty() || !bwd_queue.empty()) {
        if (is_timeout()) break;
        if (!fwd_queue.empty()) {
            auto fwd_frontier = fwd_queue.top();
            fwd_queue.pop(); // This clears the top entry.
            //fwd_processed_path[fwd_frontier.path.back()].insert(fwd_frontier);
            // Matched for this frontier.

            if (check_for_complete_path(fwd_frontier.path.back(), bwd_processed_path)) {
                // found one

                auto patt1 = getTime();
                pattern_evaluation(gm, target, fwd_frontier, bwd_processed_path);
                auto patt2 = getTime();
                pattern_evaluation_time += getInterval(patt1, patt2);

            }
            // Forward BFS expansion.
            if ((fwd_frontier.path.size() + 1) / 2 >= fwd_path_length_limit)
                continue;
            auto bfs1 = getTime();
            BFS_expansion(gm, fwd_frontier, fwd_queue, fwd_processed_path,
                          target, context, fwd_delayed_queue);
            auto bfs2 = getTime();
            BFS_expansion_time += getInterval(bfs1, bfs2);
        }

        if (!bwd_queue.empty()) {
            auto bwd_frontier = bwd_queue.top();
            bwd_queue.pop(); // This clears the top entry.

            //bwd_processed_path[bwd_frontier.path.back()].insert(bwd_frontier);

            // Matched for this frontier
            if (check_for_complete_path(bwd_frontier.path.back(), fwd_processed_path)) {
                auto patt1 = getTime();
                pattern_evaluation(gm, target, bwd_frontier, fwd_processed_path, false);
                auto patt2 = getTime();
                pattern_evaluation_time += getInterval(patt1, patt2);
            }

            // Backward BFS expansion.
            if ((bwd_frontier.path.size() + 1) / 2 >= bwd_path_length_limit)
                continue;

            auto bfs1 = getTime();
            BFS_expansion(gm, bwd_frontier, bwd_queue, bwd_processed_path,
                          target, context, bwd_delayed_queue, false);

            auto bfs2 = getTime();
            BFS_expansion_time += getInterval(bfs1, bfs2);

        }
    }
}


void
FMiner::BFS_expansion(GraphManager &gm, const BFS_entry &frontier,
                      BFS_queue &bfs_queue, node_path &processed_path,
                      uint32_t target, uint32_t context,
                      BFS_queue &delayed_queue,
                      bool is_fwd) {
    auto frontier_node = frontier.path.back();
    if (frontier.path.size() > 1 && (frontier_node == context || frontier_node == target))
        return;
    size_t edge_start = gm.nodes[frontier_node];
    size_t adj_sz = gm.nodes[frontier_node + 1] - edge_start;


    // Avoid repetition.
    unordered_set<uint32_t> nodes_on_path;
    for (int i = 0; i < frontier.path.size(); i += 2)
        nodes_on_path.insert(frontier.path[i]);

    for (int i = 0; i < adj_sz; ++i) {
        auto edge = gm.edges[edge_start + i];
        uint32_t edge_id = extractHigh32bits(edge);
        uint32_t end_node = extractLow32bits(edge);
        int end_node_adj = gm.nodes[end_node + 1] - gm.nodes[end_node];
        // Ignore nodes with no type
        if (gm.nid2types[end_node].empty()) continue;

        // check for path length_limit
        if (nodes_on_path.size() + 1 == path_length_limit) {
            if (is_fwd && end_node != context) continue;
            if (!is_fwd && end_node != target) continue;
        }

        if (nodes_on_path.count(end_node)) continue; // Cycle detected.
        if (!node_score.count(end_node))
            continue; // irrelevant

        // needs to check the k-th score in top-k queue.
        double current_path_score = 1;
        if (is_fwd) {
            // in forward direction.
            auto _exp = (double) nodes_on_path.size();
            double current_decay_factor = pow(decay_factor, -_exp);
            current_path_score = min(frontier.score,
                                     node_score[end_node] * current_decay_factor);
        } else {
            // in backward direction.
            for (int decay_exp = 0; decay_exp < frontier.path.size(); decay_exp += 2) {
                int decay_exp_reverse = (int) nodes_on_path.size() - decay_exp / 2 - 1;
                double current_decay_factor = pow(decay_factor, -decay_exp_reverse);
                current_path_score = min(current_path_score,
                                         node_score[frontier.path[decay_exp]] * current_decay_factor);
            }

            current_path_score = min(current_path_score, node_score[end_node]);
        }


        if (!above_k_th_pattern_score(current_path_score)) {
            // Should compare with the k-th highest score
            // This has a strong power when used in delayed search.
            continue;
        }

        vector<uint32_t> derived_path(frontier.path);

        if (is_fwd)
            derived_path.push_back(edge_id);
        else // reverse edge. No need bother with edge direction later on.
            derived_path.push_back(gm.reverseEdgeTypeId(edge_id));

        derived_path.push_back(end_node);
        BFS_entry tmp_entry(current_path_score, derived_path);
//        if (processed_path.count(end_node) && processed_path[end_node].count(tmp_entry)) {
//            continue; // should be impossible. Becomes impossible in the second round processing.
//        }
        /// If the degree is too high, we should delay processing it
        if (this->delay_highdegree_expansion && end_node_adj >= degree_threshold) {
            delayed_queue.emplace(tmp_entry);
        } else {
            bfs_queue.emplace(tmp_entry);
        }
        processed_path[end_node].insert(tmp_entry);
    }
}

void
FMiner::pattern_evaluation(GraphManager &gm,
                           uint32_t target,
                           const BFS_entry &entry,
                           node_path &processed_partial_path,
                           bool is_fwd) {
    for (auto &processed_entry : processed_partial_path[entry.path.back()]) {
        if (is_timeout()) break;
        // Test path length limit. Minus 1 because of one common intersect node.
        int path_len = int(processed_entry.path.size() + entry.path.size()) / 2;
        if (path_len > this->path_length_limit)
            continue;

        auto path_score1 = getTime();
        double score;
        vector<uint32_t> path;
        if (is_fwd) {
            score = min(entry.score, processed_entry.score * pow(decay_factor, -entry.path.size()));
            path = entry.path;
            for (int i = 1; i < processed_entry.path.size(); ++i) {
                int idx = (int) processed_entry.path.size() - i - 1;
                path.push_back(processed_entry.path[idx]);
            }
        } else {
            score = min(processed_entry.score, entry.score * pow(decay_factor, -processed_entry.path.size()));
            path = processed_entry.path;
            for (int i = 1; i < entry.path.size(); ++i) {
                int idx = (int) entry.path.size() - i - 1;
                path.push_back(entry.path[idx]);
            }
        }
        auto path_score2 = getTime();
        pattern_eval_path_score_time += getInterval(path_score1, path_score2);

        // Avoid repeatedly processing a path.
        if (is_path_visited(path)) {
            visited_path_num++;
            continue;
        }

        // Check the path has no cycle O(n).
        if (has_cycle(path)) continue;

        if (!above_k_th_pattern_score(score)) {
            // Should test against top-k result.
            top_k_pruned_path_num++;
            continue;
        }

        // Start finding matching node sets here.
        // This can be optimized by memorize intermediate result.
        // First, get the basic path pattern.
        auto node_type1 = getTime();
        bool skip = false;
        vector<string> path_pattern;
        for (int i = 0; i < path.size() - 1 && !is_timeout(); i += 1) {
            string most_rel_type;
            if (i % 2 == 0) {
                if (get_node_type(gm, target, path[i], most_rel_type)) {
                    path_pattern.push_back(most_rel_type);
                } else {
                    skip = true;
                    break;
                }
            } else {
                path_pattern.push_back(gm.typeId2Name[path[i]]);
            }
        }
        auto node_type2 = getTime();
        pattern_eval_node_type_time += getInterval(node_type1, node_type2);
        if (skip) continue;
        path_pattern.push_back(to_string(path.back()));

        if (gm.type2nid[path_pattern[0]].size() < num_peer_entities) {
            continue;
        }

        // Second, enumerate all path patterns
        vector<vector<string>> possible_patterns;
        auto penum1 = getTime();
        pattern_enumeration(path, path_pattern, 0, possible_patterns);
        auto penum2 = getTime();
        pattern_enumeration_time += getInterval(penum1, penum2);

        // Third, find matching node sets for each patterns and insert patterns to the top-k queue.
        bool not_significant_anymore = false;
        if (verbose)
            cout << "-----------------------------" << endl;
//        cout << "Total size = " << possible_patterns.size() << endl;
        auto spe1 = getTime();
        for (int k = 0; k < possible_patterns.size(); ++k) {
            auto &per_pattern = possible_patterns[k];
            /// Remove repetitions below

            string _pattern_key;
            pattern_memo::node_pattern_sig(path, per_pattern, path.back(),
                                      path.size() - 1, _pattern_key);

            if (processed_patterns.count(_pattern_key)) {
                continue;
            }

            processed_patterns.insert(_pattern_key);

            if (not_significant_anymore) {
                continue;
            }

            /// Remove repetitions above.

            bool size_exceed = false;
            double size_limit_time1 = size_limit_time;
            last_level_matching_time = 0;
            vector<pair<unordered_set<uint32_t>, double>> matching_nodes((path.size() + 1) / 2);
            use_size_bound = false;
            auto mnode_set1 = getTime();
            pattern_evaluation_num++;
            matching_node_set(gm, path, per_pattern, path.back(),
                              matching_nodes,
                              (int) path.size() - 1,
                              size_exceed);
            auto mnode_set2 = getTime();
            matching_node_set_call_time += getInterval(mnode_set1, mnode_set2);

            if (verbose)
                printf("######\n node set time = %.3lf ms, size limit time = %.3lf, last_level_matchingTime = %.3lf ms\n",
                       getInterval(mnode_set1, mnode_set2), size_limit_time - size_limit_time1, last_level_matching_time);

            if (verbose) {
                for (int i = 0; i < matching_nodes.size(); ++i) {
                    cout << matching_nodes[i].first.size() << " <|> ";
                }
                cout << " --\t-- ";
                for (int i = 0; i < per_pattern.size(); ++i) {
                    cout << per_pattern[i] << ", ";

                }
                string full_pkey;
                pattern_memo::pattern_sig(path, per_pattern,
                                          path.back(), path.size() - 1,
                                          full_pkey);
                cout << " --\t-- " << full_pkey;
                cout << endl;
            }

            if (size_exceed) {
                if (verbose)
                    cout << "~~~~~~~~~~~~~~~~~~~~~ size_exceed\n######" << endl;
                continue;
            }

            if (matching_nodes[0].first.size() < num_peer_entities) {
                not_enough_peers++;
                // The rest of course fails
                if (verbose)
                    printf("###### not significant. pattern_0 possible size %ld \n", gm.type2nid[per_pattern[0]].size());
                not_significant_anymore = true;

//                unusable_prefix_pattern[per_pattern[0]] = per_pattern[1];

                continue;
            }

            // Insert to top-k heap. Calculating COF is at the end of the processing.
            double pattern_relevant_score = 1;

            for (int i = 0; i < matching_nodes.size(); ++i) {
                double matching_node_score = matching_nodes[i].second / matching_nodes[i].first.size();
                matching_node_score *= pow(decay_factor, -i);
                pattern_relevant_score = min(pattern_relevant_score, matching_node_score);
            }

            // compare with path score
            pattern_relevant_score = min(pattern_relevant_score, score);
            auto insert_start = getTime();

            insert_topk_pattern(pattern_relevant_score,
                                path,
                                per_pattern,
                                matching_nodes[0].first);

            auto insert_end = getTime();
            topk_insertion_time += getInterval(insert_start, insert_end);
//            printf("\nAbove time, patternEval=%.3lf ms, ScoreCal=%.3lf ms, score=%.8lf, visited_path=%d, pruned_path=%d\n",
//                   getInterval(mnode_set1, mnode_set2), getInterval(mnode_set2, score_finished), pattern_relevant_score,
//                   visited_path_num, top_k_pruned_path_num);

        }
        auto spe2 = getTime();
        pattern_eval_time2 += getInterval(spe1, spe2);

    }


}


bool
FMiner::matching_node_set(GraphManager &gm,
                          std::vector<uint32_t> &path,
                          std::vector<std::string> &pattern,
                          uint32_t frontier,
                          std::vector<std::pair<std::unordered_set<uint32_t>, double>> &matching_nodes,
                          int idx,
                          bool &size_exceeded) {
    if (is_timeout()) {
        return false;
    }
    int matching_node_idx = idx / 2;
    if (idx == 0) {
        // Check node type here
        // Note that the node type has already been checked.
        if (!matching_nodes[matching_node_idx].first.count(frontier)) {
            matching_nodes[matching_node_idx].first.insert(frontier);
            matching_nodes[matching_node_idx].second += get_node_score(frontier);

            // impossible for context node reach here

            if (size_bound_opt_on) {
                int size_limit = 0; // no use here
                size_exceeded = exceed_size_limit(gm,
                                                  matching_node_idx,
                                                  pattern[matching_node_idx],
                                                  matching_nodes[matching_node_idx].first.size(),
                                                  path[idx + 1],
                                                  size_limit);

            }

        }

        return true;
    }

    /// USE memo for optimizing different patterns.
    string _pattern_key;
    if (use_memo) {
        pattern_memo::pattern_sig(path, pattern, frontier, idx, _pattern_key);
        if (dp_pattern_memo.has_key(_pattern_key)) {
            memo_use_num += 1;
            auto mn_sets = dp_pattern_memo.get_by_key(_pattern_key);
            if (mn_sets.empty())
                return false;

            // insert the frontier itself
            if (!matching_nodes[matching_node_idx].first.count(frontier)) {
                matching_nodes[matching_node_idx].first.insert(frontier);
                matching_nodes[matching_node_idx].second += get_node_score(frontier);
            }
            int limit_size = std::numeric_limits<int>::max();
            if (size_bound_opt_on) {
                /// memo takes effect for size bound
                if (dp_pattern_memo.is_exceed(_pattern_key)) {
                    size_exceeded = true;
//                    if (verbose)
//                        printf("-- First size_exceed level: %d\n", )
                    return true;
                }


                if (size_bound_opt_on) {
                    size_exceeded = exceed_size_limit(gm, matching_node_idx,
                                                      pattern[matching_node_idx],
                                                      matching_nodes[matching_node_idx].first.size(),
                                                      path[matching_node_idx + 1],
                                                      limit_size);
                    if (size_exceeded) return true;
                }
            }

            // reuse nodes before this node
//            bool tmp_size_exceed = size_exceeded;
            for (int i = 0; i < mn_sets.size(); ++i) {
                for (auto n : mn_sets[i]) {
                    if (matching_nodes[i].first.size() > limit_size) {
                        size_exceeded = true;
                        break;
                    }
                    if (!matching_nodes[i].first.count(n)) {
                        matching_nodes[i].first.insert(n);
                        matching_nodes[i].second += get_node_score(n);
                    }
                }
            }
//            size_exceeded = tmp_size_exceed;

            return true;
        }
    }

    uint32_t edge_type = path[idx - 1];
    uint32_t reversed_edge_type = gm.reverseEdgeTypeId(edge_type);

    // handle next node type being a node
    string next_node_type = pattern[matching_node_idx - 1];
    int instance_node = node_type_id(next_node_type);
    if (instance_node != -1) {
        bool matched = matching_node_set(gm, path, pattern, instance_node,
                                         matching_nodes, idx - 2,
                                         size_exceeded);

        if (matched) {
            matching_nodes[matching_node_idx].first.insert(frontier);
            matching_nodes[matching_node_idx].second += get_node_score(frontier);
            // suppose This frontier corresponds to node type human.
            // next_instance_node<-human<- sexOfHuman
            // There is no check if don't do the following.
            if (size_bound_opt_on) {
                int size_limit = 0; // no use here
                size_exceeded = exceed_size_limit(gm,
                                                  matching_node_idx,
                                                  pattern[matching_node_idx],
                                                  matching_nodes[matching_node_idx].first.size(),
                                                  path[idx + 1],
                                                  size_limit);
            }
            // Test if size_exceed // no need for single node.
            return true;
        } else {
            return false;
        }
    }

    int edge_pos = gm.firstEdgePos(frontier, reversed_edge_type);
    if (edge_pos == -1) return false;

    bool has_match = false;
    /// never include the current frontier.
    vector<pair<unordered_set<uint32_t>, double>> tmp_matching_nodes(matching_node_idx);
    for (int i = edge_pos; i < (int) gm.nodes[frontier + 1] && !is_timeout(); ++i) {
        uint32_t eid = extractHigh32bits(gm.edges[i]);

        if (eid != reversed_edge_type) break;
        uint32_t end_node = extractLow32bits(gm.edges[i]);
        if (matching_node_idx >= 2) { // look ahead
            uint32_t next_edge_type = path[idx - 3];
            uint32_t reversed_next_edge_type = gm.reverseEdgeTypeId(next_edge_type);
            if (!gm.typeId2Count[reversed_next_edge_type].count(end_node)) {
                continue;
            }
        }
        // Should not do this (otherwise the pattern matching instances are incomplete):
        // if (!node_score.count(end_node)) continue;

        // Internal node cannot be context nodes.
//        if (context_nodes.count(end_node)) continue;

        // End node does not match the required type
        if (!node_has_type(gm, end_node, next_node_type))
            continue;

        // end node has already been processed
        if (!use_memo && matching_nodes[matching_node_idx - 1].first.count(end_node))
            continue;

        if (matching_node_set(gm, path, pattern, end_node,
                              tmp_matching_nodes, idx - 2,size_exceeded)) {
            has_match = true;
            // test size bound. If no match, then size_exceeded cannot be true.
            if (size_bound_opt_on) {
                if(size_exceeded) break; // do not quit so fast. We need to fill memo, so to avoid future repetition.
            }
        }
    }

    if (has_match) {
        if (!matching_nodes[matching_node_idx].first.count(frontier)) {
            matching_nodes[matching_node_idx].first.insert(frontier);
            matching_nodes[matching_node_idx].second += get_node_score(frontier);
        }
        if (size_bound_opt_on) {
            int limit_size_current_node = std::numeric_limits<int>::max();
            if (!size_exceeded) { // avoid overwrite the size_exceed in case it is already true.
                size_exceeded = exceed_size_limit(gm, matching_node_idx, pattern[matching_node_idx],
                                                  matching_nodes[matching_node_idx].first.size(),
                                                  path[matching_node_idx + 1],
                                                  limit_size_current_node);
            }
        }

        for (int i = 0; i < tmp_matching_nodes.size(); ++i) {
            if (size_bound_opt_on && size_exceeded) break;
            int limit_size = std::numeric_limits<int>::max();
            if (size_bound_opt_on) {
                size_exceeded = exceed_size_limit(gm, i, pattern[i],
                                                  tmp_matching_nodes[i].first.size(),
                                                  path[i + 1],
                                                  limit_size);
                if (size_exceeded) {
                    // we need to fill memo. Do not return so fast.
                    break;
                }
            }

            // merge and derive new memo
            for (auto n : tmp_matching_nodes[i].first) {
                if (!matching_nodes[i].first.count(n)) {
                    matching_nodes[i].first.insert(n);
                    matching_nodes[i].second += get_node_score(n);
                    if (size_bound_opt_on && matching_nodes[i].first.size() > limit_size) {
                        size_exceeded = true;
                        break;
                    }
                }

            }
        }

        if (use_memo) {
            // save it to memo
            if (size_bound_opt_on)
                dp_pattern_memo.insert(_pattern_key, tmp_matching_nodes, size_exceeded);
            else
                dp_pattern_memo.insert(_pattern_key, tmp_matching_nodes);
        }
        return true;

    } else {
        // insert empty
        if (use_memo) {
            dp_pattern_memo.insert(_pattern_key);
        }

        return false;
    }
}

/// most relevant node type: used for avoiding repetitive computation
bool
FMiner::get_node_type(GraphManager &gm,
                      uint32_t target,
                      uint32_t node,
                      std::string &most_rel_type) {

    if (most_rel_node_type.count(node)) {
        // Already calculated.
        if (most_rel_node_type[node].empty()) {
            return false;
        }
        most_rel_type = most_rel_node_type[node];
        return true;
    }

    const auto &nodeTypes = gm.nid2types[node];
    if (nodeTypes.empty()) {
        return false;
    }

    // impossible for score to be larger than 1
    double max_type_rel_score = 0;
    for (const auto &t : nodeTypes) {
        double tmp_score = 0;
        if (type_rel_score.count(t)) {
            tmp_score = type_rel_score[t];
        } else {
            // consider this node type
            double accum_count = 0;
            for (uint32_t nid : gm.type2nid[t]) {
                if (nid != target) {
                    accum_count += 1;
                    if (node_score.count(nid)) {
                        tmp_score += node_score[nid];
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
        most_rel_node_type[node] = "";
        return false;
    }
    most_rel_node_type[node] = most_rel_type;
    return true; // process okay
}


void
FMiner::pattern_enumeration(vector<uint32_t> &path,
                            vector<std::string> &path_pattern,
                            size_t idx,
                            vector<vector<std::string>> &enumerated_path_pattern) {
//    vector<string> tmp_pattern;
//    for (int i = 0; i < path_pattern.size(); i+=2) {
//        tmp_pattern.emplace_back(path_pattern[i]);
//    }
//    enumerated_path_pattern.emplace_back(tmp_pattern);
//
//    for (int i = (int)tmp_pattern.size() - 2; i > 0 ; --i) {
//        enumerated_path_pattern.emplace_back(tmp_pattern);
//        enumerated_path_pattern.back()[i] = to_string(path[i * 2]);
//    }


    // With cut-off node, we do not need the following.
    if (idx == path_pattern.size() - 1) {
        enumerated_path_pattern.back().emplace_back(path_pattern.back());
        return;
    }
    if (idx == 0) {
        enumerated_path_pattern.emplace_back();
        enumerated_path_pattern.back().emplace_back(path_pattern[0]);
        pattern_enumeration(path, path_pattern, idx + 2, enumerated_path_pattern);
    } else {
        auto tmp = enumerated_path_pattern.back();
        enumerated_path_pattern.back().emplace_back(path_pattern[idx]);
        pattern_enumeration(path, path_pattern, idx + 2, enumerated_path_pattern);
        tmp.emplace_back(to_string(path[idx]));
        enumerated_path_pattern.emplace_back(tmp);
        pattern_enumeration(path, path_pattern, idx + 2, enumerated_path_pattern);
    }
}


void
FMiner::topk_cof_extraction(GraphManager &gm,
                            uint32_t target) {
    // Should avoid processing the same set of peers.
    final_topk_pattern_entries.resize(topk_patterns.size());
    size_t idx_place = final_topk_pattern_entries.size() - 1;
    while (!topk_patterns.empty()) {
        size_t peer_size = topk_patterns.top().peer_entities.size();
        auto tmp_entry = topk_patterns.top();
//        if (visited_peer_set.count(peer_size)) {
//            final_topk_pattern_entries[idx_place--] = std::move(tmp_entry);
//            topk_patterns.pop();
//            continue;
//        }
        visited_peer_set.insert(peer_size);
        cof_scoring(gm, target, tmp_entry);
        final_topk_pattern_entries[idx_place--] = std::move(tmp_entry);
        topk_patterns.pop();
    }
}


void
FMiner::cof_scoring(GraphManager &gm, uint32_t target, topk_entry &entry) {
    auto &peer_entities = entry.peer_entities;
    int edge_pos = gm.nodes[target];
    int adj_sz = gm.nodes[target + 1] - edge_pos;

    unordered_map<uint32_t, bool> legal_edge_ids;
    topk_cof_queue top_k_cofs;

    int reason1 = 0, reason2 = 0, reason3 = 0;
    /// Important (affecting result quality): What attributes are interesting enough for a probe?
    /// Some attributes are striking but not interesting enough.
    for (int i = 0; i < adj_sz; ++i) {
        uint64_t edge = gm.edges[edge_pos + i];
        uint32_t edge_tid = extractHigh32bits(edge);
        uint32_t end_node = extractLow32bits(edge);

        /// node relevance
//        if (!node_score.count(end_node)) {
//            reason1++;
//            continue;
//        }

        /// node importance
//        uint32_t degree_end_node = gm.nodes[end_node + 1] - gm.nodes[end_node];
//        if (degree_end_node < degree_importance_limit) {
//            reason2++;
//            continue;
//        }

        // edge eligibility
        if (!legal_edge_ids.count(edge_tid))
            legal_edge_ids[edge_tid] = validate_edge(gm, peer_entities, edge_tid, end_node);
        if (!legal_edge_ids[edge_tid]) {
            reason3++;
            continue;
        }

        // calculate strikingness score
        unordered_map<uint32_t, unordered_set<uint32_t>> attribute_peers;
        unordered_map<uint32_t, unordered_set<uint32_t>> peer_attributes;
        //unordered_map<uint32_t, double> attribute_counts; // used for calculating striking score
        double end_node_freq = 1;
        unordered_map<uint32_t, double> val_frequency;
        val_frequency[end_node] = 1;

        for (auto p : peer_entities) {
            if (p == target) continue;
            auto _key = combine2u32(p, edge_tid);
            if (!gm.snidEid2pos.count(_key)) continue;
            auto peer_edge_pos = gm.snidEid2pos[_key];
            if (peer_edge_pos == -1) continue;

//            auto expected_hit = combine2u32(edge_tid, end_node);
//            bool has_end_node = std::binary_search(gm.edges + peer_edge_pos, gm.edges + gm.nodes[p+1], expected_hit);
//            if (has_end_node) {
//                val_frequency[end_node] += 1;
//            }
            bool has_end_node = false;
            unordered_set<uint32_t> tmp_attrs;
            uint32_t peer_val;
            for (int j = peer_edge_pos; j < gm.nodes[p + 1] && !has_end_node; ++j) {
                auto attribute_val = extractLow32bits(gm.edges[j]);
                auto attribute_edge = extractHigh32bits(gm.edges[j]);
                if (attribute_edge != edge_tid) break;
                if (attribute_val == end_node) {
                    has_end_node = true;
                    peer_val = end_node;
                    break;
                }


                if (j == peer_edge_pos) {
                    peer_val = attribute_val;
                }
                // Avoid deletion from the map if has end_node.
                //attribute_peers[attribute_val].insert(p);
                tmp_attrs.insert(attribute_val);
                peer_attributes[p].insert(attribute_val);
            }
//            cout << "---------" << endl;

            if (val_frequency.count(peer_val)) {
                val_frequency[peer_val] += 1;
            } else {
                val_frequency[peer_val] = 1;
            }

            if (!has_end_node) {
                for (auto attr : tmp_attrs) {
                    attribute_peers[attr].insert(p);
                }
            }
        }

        std::vector<std::pair<uint32_t, int>> tmp_exp;
        double higher_than_end_node_freq = 0;
        int total_different_attr = 1; // +1 for the target attribute end_node
        double total_freq = 0;
        double highest_freq = 0;
        uint32_t highest_node = end_node;
        end_node_freq = val_frequency[end_node];
        for (auto &kv : val_frequency) {
            if (kv.second > end_node_freq) {
                higher_than_end_node_freq += kv.second;
            }
            if (kv.second > highest_freq) {
                highest_freq = kv.second;
                highest_node = kv.first;
            }
            total_freq += kv.second;
        }


        // Considering from most frequent attributes to least ones.
        // Will it be more interesting to consider the opposite way.
//        while (true && !is_timeout()) {
//            bool updated = false;
//            int highest_attr, highest_peer_num = 0;
//            for (auto &attr_peers : attribute_peers) {
//                if (attr_peers.second.size() > highest_peer_num) {
//                    highest_peer_num = (int) attr_peers.second.size();
//                    highest_attr = (int) attr_peers.first;
//                    updated = true;
//                }
//            }
//            if (!updated) break;
//            if (highest_peer_num >= end_node_freq) {
////                tmp_exp.emplace_back(make_pair(highest_attr, highest_peer_num));
//                tmp_exp.emplace_back(highest_attr, highest_peer_num);
//                higher_than_end_node_freq += highest_peer_num;
//            }
//            total_different_attr++;
//
//            for (auto to_del_peer: attribute_peers[highest_attr]) {
//                for (auto modify_attr: peer_attributes[to_del_peer]) {
//                    if (modify_attr != highest_attr) {
//                        attribute_peers[modify_attr].erase(to_del_peer);
//                    }
//                }
//            }
//            attribute_peers[highest_attr].clear();
//        }
        if (is_timeout())
            return;

        // The distribution is in attribute_counts.
        // get COF score. The larger the better.
        //double current_cof_score = higher_than_end_node_freq / peer_entities.size();
        double current_cof_score = higher_than_end_node_freq / total_freq;
        cof_entry _tmp_cof(current_cof_score, (int) end_node_freq, edge_tid, end_node, tmp_exp);

        _tmp_cof.highest_freq = highest_freq;
        _tmp_cof.highest_IRI = gm.nid2IRI[highest_node];

        if (top_k_cofs.size() < cof_topk_num) {
            top_k_cofs.emplace(_tmp_cof);
        } else {
            if (top_k_cofs.top().striking_score < current_cof_score) {
                top_k_cofs.pop();
                top_k_cofs.emplace(_tmp_cof);
            }
        }
    }
    // printf("Irrelevant: %d, low degree: %d, illegal edge: %d\n", reason1, reason2, reason3);
    // Storing top-k COFs.
    entry.top_k_cofs.resize(top_k_cofs.size());
    size_t idx_place = top_k_cofs.size() - 1;
    while (!top_k_cofs.empty()) {
        entry.top_k_cofs[idx_place--] = top_k_cofs.top();
        top_k_cofs.pop();
    }
}


bool
FMiner::validate_edge(GraphManager &gm, std::unordered_set<uint32_t> &peers, uint32_t edge_tid, uint32_t end_node) {
    /// First, guarantee the edge id appears in most peer entities.
    int count = 0;
    for (auto &p : peers)
        count += gm.typeId2Count[edge_tid].count(p);

    if (count < attribute_occurrence_ratio_limit * peers.size()) {
        // printf("Failed: Have: %ld, Got: %ld.\n", peers.size(), count);
        return false;
    }
//    printf("Succeed: Have: %ld, Got: %d.\n", peers.size(), count);
    return true;
}
