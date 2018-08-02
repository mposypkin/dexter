/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   bnbseq.hpp
 * Author: mposypkin
 *
 * Created on July 31, 2018, 5:36 PM
 */

#ifndef BNBSEQ_HPP
#define BNBSEQ_HPP

#include <iostream>
#include <vector>
#include <expression/expr.hpp>


#include "bnbcore.hpp"

/**
 * Performs iterations
 * @param es equality system
 * @param boxlist the box list
 * @param mind minimum diameter of a box
 * @param maxSteps the maximal number of steps to perform
 */
void iterate(const EqualitySystem& es, std::vector<Box>& boxlist, std::vector<Box>& boundary, std::vector<Box>& internal, double mind, long long int &maxSteps) {
    while (!boxlist.empty()) {
        Box b = boxlist.back();
        boxlist.pop_back();
        BoxSort sort = sortBox(es, b);
        switch (sort) {
            case BoxSort::boundary:
                if (maxEdge(b) <= mind) {
                    boundary.push_back(b);
                } else {
                    splitBox(b, boxlist);
                }
                break;
            case BoxSort::external:
                break;
            case BoxSort::internal:
                internal.push_back(b);
        }
    }
}




#endif /* BNBSEQ_HPP */

