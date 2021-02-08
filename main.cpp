#include <iostream>
#include "GraphManager.h"
#include "RandSim.h"
#include "ExceptionalFactMining.h"
#include "settings.h"
#include "QueryInterface.h"


using namespace std;

void find_cofs(int h_num = 2,
               int top_k_pattern = 20,
               int length_limit = 4,
               bool reverse_search = true,
               bool optDFS = true);

void testSingleQuery(GraphManager &gm);

void testBlindSearch(GraphManager &gm);



int main(int argc, char **argv) {
    std::cout << "Hello, yueji!" << std::endl;


    //Search
    int h_num = 2;
    if (argc == 2)
        h_num = stoi(argv[1]);

    // default
    int top_k_pattern = 20;
    int length_limit = 4;
    bool reverse_search = true;
    bool optDFS = true;

    if (argc == 1 + 4) {
        // three args:
        // top-k pattern
        top_k_pattern = stoi(argv[1]);
        // path length limit
        length_limit = stoi(argv[2]);
        // 1 for opt reversed search, 0 for naive search
        reverse_search = (bool) stoi(argv[3]);
        // 1 for opt pattern eval, 0 for naiveDFS
        optDFS = (bool) stoi(argv[4]);
    }


    find_cofs(h_num, top_k_pattern,
              length_limit,
              reverse_search,
              optDFS);

    printf("End of Processing.\n");
    return 0;
}


void testBlindSearch(GraphManager &gm, int target) {
    ExceptionalFactMining efm;
    efm.global_topk = 50;
    efm.enumPathPatternBlind(gm, target);

    exit(123);
}


void testSingleQuery(GraphManager &gm,
                     int top_k_pattern = 20,
                     int length_limit = 4,
                     bool reverse_search = true,
                     bool optDFS = true) {

    QueryInterface qi;

    /// Find nodeId from IRI
    /// Kamala Harris (Q10853588) = 891030, Mike Pence (Q24313) = 2393228
//    uint32_t harris = 891030;
    //    qi.getQueryNidPairs(gm, "../../testData/Harris.txt");
//    for (auto &p : qi.queryNidPairs) {
//        cout << p.first << "," << p.second << endl;
//    }
//    testBlindSearch(gm, harris);
    //exit(124);

#ifdef ServerRun
    // below is used for real data
//    qi.processAllQueries(gm, "../../testData/first2000Pairs.txt",
//                         top_k_pattern,
//                         length_limit,
//                         reverse_search,
//                         optDFS);

    /// For new FMiner algorithm
//    qi.process_all_queries_fminer(gm, "../../testData/first2000Pairs.txt",
//    qi.process_all_queries_fminer(gm, "../../testData/testHard.txt",
//                                  top_k_pattern,
//                                  length_limit,
//                                  reverse_search,
//                                  optDFS);

    // delay_highdegree_expansion+, use_memo+, size_bound_opt_on+
    qi.process_all_queries_fminer(gm, "../../testData/testHard.txt",
                                  top_k_pattern,
                                  length_limit,
                                  true,
                                  true,
                                  true);
    exit(124);
    // delay_highdegree_expansion+, use_memo-, size_bound_opt_on-
    qi.process_all_queries_fminer(gm, "../../testData/testHard.txt",
                                  top_k_pattern,
                                  length_limit,
                                  true,
                                  false,
                                  false);

    // delay_highdegree_expansion-, use_memo+, size_bound_opt_on-
    qi.process_all_queries_fminer(gm, "../../testData/testHard.txt",
                                  top_k_pattern,
                                  length_limit,
                                  false,
                                  true,
                                  false);

    // delay_highdegree_expansion-, use_memo-, size_bound_opt_on+
    qi.process_all_queries_fminer(gm, "../../testData/testHard.txt",
                                  top_k_pattern,
                                  length_limit,
                                  false,
                                  false,
                                  true);

    // delay_highdegree_expansion-, use_memo-, size_bound_opt_on-
    qi.process_all_queries_fminer(gm, "../../testData/testHard.txt",
                                  top_k_pattern,
                                  length_limit,
                                  false,
                                  false,
                                  false);

#else
//    qi.processAllQueries(gm, "../../testData/querySinglePair.txt",
//                         top_k_pattern,
//                         length_limit,
//                         reverse_search,
//                         optDFS); //for development

    /// For new FMiner algorithm
    qi.process_all_queries_fminer(gm, "../../testData/querySinglePair.txt",
                                  top_k_pattern,
                                  length_limit,
                                  reverse_search,
                                  optDFS);

//    qi.processAllQueries(gm, "../../testData/allPairs.txt");
#endif


//    qi.processAllQueries(gm, "../../testData/first100Pair.txt");
    exit(123);
}


void find_cofs(int h_num, int top_k_pattern,
               int length_limit,
               bool reverse_search,
               bool optDFS) {
    GraphManager gm;

#ifdef  ServerRun
    //trump is 276748, clinton is 276746
    string edgeFile = "/home/yangyueji/Data/graphsWithoutInstanceOf/edgesWithoutInstanceOf.csv";
    string node2typeFile = "/home/yangyueji/Data/graphsWithoutInstanceOf/nodeTypesNames.csv";
    string serializeFile = "/home/yangyueji/Data/graphsWithoutInstanceOf/graphAllEdges.dat";
#else
    string edgeFile = "../../testData/sampledGraphWOInstOf/sampledEdges.csv";
    string node2typeFile = "../../testData/sampledGraphWOInstOf/sampledNodeTypes.csv";
    string serializeFile = "../../testData/sampledGraphWOInstOf/graph.dat";
#endif


#ifdef Deserialize
    gm.deserializeFromDisk(serializeFile.c_str());
#else
    gm.readEdges(edgeFile.c_str());
    gm.readNodeTypes(node2typeFile.c_str());
    gm.setEdgeType2Pos();
    gm.setIRI2NidMap();
    gm.sortSameTypeNodesByDegree();
    gm.serializeToDisk(serializeFile.c_str());
    exit(123);
#endif


    printf("Node Num = %ld, edge num = %ld.\n", gm.node_sz, gm.edge_sz);
    // find node ids given IRI
    QueryInterface qi;
    string nidQueryFile = "../../testData/userStudy.txt";
    qi.processAllQueries(gm, nidQueryFile);
    exit(125);

    qi.process_all_queries_fminer(gm, "../../testData/userStudy.txt",
                                  top_k_pattern,
                                  length_limit,
                                  true,
                                  true,
                                  true);

    exit(125);

    testSingleQuery(gm,
                    top_k_pattern,
                    length_limit,
                    reverse_search,
                    optDFS);
//    testMultiQueries(gm);



// Hillary Clinton Q6294, nid = 3198,
// Trump Q22686, nid = 3199

///find Hillary and Trump
//    for (int i = 0; i < gm.node_sz; ++i) {
//        auto s = gm.nid2IRI[i];
//        if (s.substr(s.length()-5) == "Q6294" || s.substr(s.length() - 6) == "Q22686") {
//            cout << i << "," << s << endl;
//        }
//    }
    ExceptionalFactMining efm;
    std::unordered_map<uint32_t, double> source2weight;
    std::unordered_map<uint32_t, double> node2score;

#ifdef ServerRun
    //    efm.contextNodes.insert(276746); // Hillary
    //    source2weight[276746] = 1;

    efm.contextNodes.insert(39);//USA
    source2weight[39] = 1;

    //    source2weight[276748] = 1;

    uint32_t entityInterest = 276748;
    efm.num_peerEntity = 20;
    efm.global_topk = 20;
    efm.degree_threshold = 10000;
    double alpha = 0.1;
    //    efm.end_node_degree_limit = 100;
#else
    //    source2weight[3199] = 1;
    //    source2weight[172637] = 1; //USA
    //    efm.contextNodes.insert(172637);

    //    source2weight[27] = 1; //new york city
    //    efm.contextNodes.insert(27);

    source2weight[3198] = 1; //Hillary
    efm.contextNodes.insert(3198);

    uint32_t entityInterest = 3199;
    //    efm.contextNodes.insert(3199);
    efm.num_peerEntity = 0;
    efm.degree_threshold = 1000;
    double alpha = 0.01;
    //    efm.end_node_degree_limit = 5;
#endif
    testBlindSearch(gm, entityInterest);

    exit(123);
}
