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

#include <cstdlib>
#include <iostream>
#include "dexter.hpp"
#include "bnbseq.hpp"



/*
 * 
 */
int main(int argc, char** argv) {
    DexterKinematicEquations de;
   
    std::vector<Interval<double>> var = { {0, 1 }, {0, 1 }, {0, 1 }, {0, 1 }, {0, 1 }, {0, 1 } };
    
    std::cout << calcInterval(de.mIG[0], var ) << "\n";
    return 0;
}

