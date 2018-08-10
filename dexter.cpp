/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   dexter.cpp
 * Author: mposypkin
 *
 * Created on July 31, 2018, 3:40 PM
 */

#include <math.h>
#include <cstdlib>
#include <iostream>

#include "dexter.hpp"
#include "bnbseq.hpp"
#include "bnbmultiset.hpp"
#include "bnbmerge.hpp"

/*
 * 
 */
int main(int argc, char** argv) {
    DexterKinematicEquations de;

    std::vector<Interval<double>> var = {
        {-5, 5},
        {-5, 5},
        {-M_PI, M_PI},
        {-M_PI, M_PI},
        {-M_PI, M_PI},
        {-M_PI, M_PI}
    };

    //    std::cout << calcInterval(de.mIG[0], var) << "\n";
    //    std::cout << sortBox(de, var);

    std::vector<Box> boxlist;
    boxlist.push_back(var);
    std::vector<Box> boundary, internal;

    constexpr double mind = 0.1;
    long long int maxSteps = 1000000;
#if 1    
    std::vector<int> coorcomp = {0, 1};
    iterateMerge(de, coorcomp, boxlist, boundary, internal, mind, maxSteps);

    for (auto b : boundary) {
        std::cout << b << "\n";
    }

    std::cout << "# Total " << boundary.size() << " boxes\n";
#else    
    std::vector<int> coorcomp = {0, 1};
    BoxCompare bc(coorcomp);
    BoxSet bsbound(bc);
    iterateSet(de, boxlist, bsbound, mind, maxSteps);
    for (auto b : bsbound) {
        std::cout << b << "\n";
    }
    std::cout << "# Total " << bsbound.size() << " boxes\n";
#endif    
    return 0;
}

