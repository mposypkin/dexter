/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   spheres.cpp
 * Author: mposypkin
 *
 * Created on August 2, 2018, 4:15 PM
 */

#include <math.h>
#include <cstdlib>
#include <iostream>

#include "spheres.hpp"
#include "bnbseq.hpp"
#include "bnbmultiset.hpp"
#include "bnbmerge.hpp"

/*
 * 
 */
int main(int argc, char** argv) {
    Spheres ce;

    std::vector<Interval<double>> var = {
        {-10, 10},
        {-10, 10},
        {-10, 10}
    };

    //    std::cout << calcInterval(de.mIG[0], var) << "\n";
    //    std::cout << sortBox(de, var);

    std::vector<Box> boxlist;
    boxlist.push_back(var);
    std::vector<Box> boundary;

    //    constexpr double mind = std::numeric_limits<double>::max();
    constexpr double mind = 0.2;
    constexpr double maxc = 4;
    long long int maxSteps = 1000000;
#if 1   
    std::vector<int> coorcomp = {0, 1};
//    const bool doLocSearch = true;
    const bool doLocSearch = false;
    iterateMerge(ce, coorcomp, boxlist, boundary, mind, maxSteps, doLocSearch);

    for (auto b : boundary) {
        std::cout << b << "\n";
    }

    std::cout << "# Total " << boundary.size() << " boxes\n";
#else    
    std::vector<int> coorcomp = {0, 1};
    BoxCompare bc(coorcomp);
    BoxSet bsbound(bc);
    iterateSet(ce, boxlist, bsbound, mind, maxSteps);
    for (auto b : bsbound) {
        std::cout << b << "\n";
    }
    std::cout << "# Total " << bsbound.size() << " boxes\n";
#endif    
    return 0;
}

