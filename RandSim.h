//
// Created by yangyueji on 8/29/19.
//

#ifndef PATHPATTERNMINING_RANDSIM_H
#define PATHPATTERNMINING_RANDSIM_H


#include "GraphManager.h"

class RandSim {
public:
    void
    randomWalkGeoSampling(const GraphManager &gm,
                          std::unordered_map<uint32_t, double> &source2weight,
                          std::unordered_map<uint32_t, double> &node2score,
                          double teleprob = 0.5);
};


#endif //PATHPATTERNMINING_RANDSIM_H
