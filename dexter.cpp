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
#include <limits>

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
        {-13, 13},
        {-16, 16},
        {-M_PI, M_PI},
        {-M_PI, M_PI},
        {-M_PI, M_PI},
        {-M_PI, M_PI}
    };

    //    std::cout << calcInterval(de.mIG[0], var) << "\n";
    //    std::cout << sortBox(de, var);

    std::vector<Box> boxlist;
    boxlist.push_back(var);
    std::vector<Box> boundary;

    constexpr double mind = 0.2;
//    constexpr double mind = 5e-1;
    long long int maxSteps = 1e8;
#if 1    
    std::vector<int> coorcomp = {0, 1};
    const bool doLocSearch = true;
//    const bool doLocSearch = false;
    iterateMerge(de, coorcomp, boxlist, boundary, mind,  maxSteps, doLocSearch);

    for (auto b : boundary) {
        std::cout << b << "\n";
    }

    std::cout << "# Total " << boundary.size() << " boxes\n";
    std::cout << "# Total " << maxSteps << " performed\n";
    
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

