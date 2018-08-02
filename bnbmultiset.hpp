/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   bnbmultiset.hpp
 * Author: mposypkin
 *
 * Created on August 2, 2018, 1:22 PM
 */

#ifndef BNBMULTISET_HPP
#define BNBMULTISET_HPP

#include <iostream>
#include <vector>
#include <expression/expr.hpp>
#include <set>

#include "bnbcore.hpp"

struct BoxCompare {

    /**
     * Constructor
     * @param compcoor vector of coordinates to compare
     */
    BoxCompare(const std::vector<int>& compcoor) : mCompCoor(compcoor) {

    }

    bool operator()(const Box &lhs, const Box &rhs) const {
        for (auto i : mCompCoor) {
            if ((lhs[i].lb() < rhs[i].lb()) || (lhs[i].rb() > rhs[i].rb()))
                return false;
        }
        return true;
    }

    std::vector<int> mCompCoor;
};


using BoxSet = std::multiset<Box, BoxCompare>;

/**
 * Performs iterations
 * @param es equality system
 * @param boxlist the box list
 * @param mind minimum diameter of a box
 * @param maxSteps the maximal number of steps to perform
 */
void iterateSet(const EqualitySystem& es, std::vector<Box>& boxlist, BoxSet& boundary, double mind, long long int &maxSteps) {
    while (!boxlist.empty()) {
        Box b = boxlist.back();
        boxlist.pop_back();
        if (boundary.upper_bound(b) == boundary.end()) {
            BoxSort sort = sortBox(es, b);
            switch (sort) {
                case BoxSort::boundary:
                    if (maxEdge(b) <= mind) {
                        boundary.insert(b);
                    } else {
                        splitBox(b, boxlist);
                    }
                    break;
                case BoxSort::external:
                    break;
                default:
                    std::cerr << "Erroneous rectangel type\n";
                    exit(-1);
            }
        }
    }
}


#endif /* BNBMULTISET_HPP */

