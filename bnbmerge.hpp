/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   bnbmerge.hpp
 * Author: mposypkin
 *
 * Created on August 8, 2018, 12:29 PM
 */

#ifndef BNBMERGE_HPP
#define BNBMERGE_HPP

#include <iostream>
#include <vector>
#include <expression/expr.hpp>


#include "bnbcore.hpp"

bool operator<=(const Box& b1, const Box& b2) {
    bool isless = true;
    const int n = b1.size();
    for (int i = 0; i < n; i++) {
        if ((b1[i].lb() < b2[i].lb()) || (b1[i].rb() > b2[i].rb())) {
            isless = false;
            break;
        }
    }
    return isless;
}

/**
 * Projects a box to a set of coordinates
 * @param pcoords coordinates for projection
 * @param box source box
 * @return new projected box
 */
Box projectBox(const std::vector<int>& pcoords, const Box& box) {
    Box nbox;
    for (auto i : pcoords) {
        nbox.push_back(box[i]);
    }
    return nbox;
}

/**
 * Evaluates merge dimension
 * @param b1 first box
 * @param b2 second box
 * @return merge dim or -1 if no such found
 */
int mergeDim(const Box& b1, const Box& b2) {
    int md = -1;
    int mdc = 0;
    const int n = b1.size();
    for (int i = 0; i < n; i++) {
        if ((b1[i].lb() == b2[i].lb()) && (b1[i].rb() == b2[i].rb())) {
            continue;
        } else if ((b1[i].rb() >= b2[i].lb()) && (b2[i].rb() >= b1[i].lb())) {
            if (mdc > 0) {
                md = -1;
                break;
            } else {
                md = i;
                mdc++;
            }
        } else {
            md = -1;
            break;
        }
    }
    return md;
}

/**
 * Merges two boxes amenable to merge
 * @param b1 first box
 * @param b2 second box
 * @return resulting box
 */
Box merge(const Box& b1, const Box& b2) {
    Box nbox;
    const int n = b1.size();
    for (int i = 0; i < n; i++) {
        double lb, rb;
        lb = std::min(b1[i].lb(), b2[i].lb());
        rb = std::max(b1[i].rb(), b2[i].rb());
        nbox.push_back(Interval<double>(lb, rb));
    }
    return nbox;
}

/**
 * Performs iterations of the B&B solver for constructing and approximation of 
 * the equality system's solution set and projecting it to the set of coordinates
 * @param es equality system
 * @param pcoords coordinates for the projection
 * @param boxlist the box list
 * @param the paving constructed
 * @param mind minimum diameter of a box
 * @param maxSteps the maximal number of steps to perform, actual number on exit
 * @param doLocSearch true if we check for root, false otherwise
 * @param maxc maximal deviation of a constraint from zero for boundary box
 */
void iterateMerge(const EqualitySystem& es, const std::vector<int>& pcoords, std::vector<Box>& boxlist, std::vector<Box>& paving, double mind, long long int &maxSteps, bool doLocSearch = true, double maxc = std::numeric_limits<double>::max()) {
    long long int counter = 0;
    while (!boxlist.empty() && (counter < maxSteps)) {
        counter ++;
        Box b = boxlist.back();
        boxlist.pop_back();
        Box bm = projectBox(pcoords, b);
        bool sameBox = false;
        std::vector<Box>::const_iterator bi;
        for (bi = paving.begin(); bi != paving.end(); bi++) {
            if (bm <= *bi) {
                sameBox = true;
                break;
            }
        }
        if (sameBox) {
            continue;
        }
        BoxSort sort = sortBox(es, b);
        switch (sort) {
            case BoxSort::boundary:
                if ((maxEdge(b) <= mind) && (getAbs(es, b) <= maxc) && ((!doLocSearch) || containsRoot(es, b))) {                    
                    for (bi = paving.begin(); bi != paving.end();) {
                        if (*bi <= bm) {
                            bi = paving.erase(bi);
                        } else if (mergeDim(*bi, bm) != -1) {
                            //std::cout << "merge " << *bi << " and " << bm << "\n";
                            bm = merge(*bi, bm);
                            paving.erase(bi);
                            bi = paving.begin();
                            continue;
                        } else {
                            bi++;
                        }
                    }
                    paving.push_back(bm);
                } else {
                    splitBox(b, boxlist);
                }
                break;
            case BoxSort::external:
                break;
            case BoxSort::internal:
                std::cerr << "Unsupported";
                exit(-1);
                break;
        }
    }
    maxSteps = counter;
}


#endif /* BNBMERGE_HPP */

