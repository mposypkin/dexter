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

/*
 * 
 */
int main(int argc, char** argv) {
    DexterKinematicEquations de;

    std::vector<Interval<double>> var = {
        {-3, 3},
        {-3, 3},
        {-M_PI, M_PI},
        {-M_PI, M_PI},
        {-M_PI, M_PI},
        {-M_PI, M_PI}
    };

    std::cout << calcInterval(de.mIG[0], var) << "\n";
    std::cout << sortBox(de, var);

    std::vector<Box> boxlist;
    boxlist.push_back(var);
    std::vector<Box> boundary, internal;

    constexpr double mind = 0.25;
    long long int maxSteps = 1000000;
    iterate(de, boxlist, boundary, internal, mind, maxSteps);

    for(auto b : boundary) {
        std::cout << b << "\n";        
    }
    
    std::cout << "Total " << boundary.size() << " boxes\n";
    return 0;
}

