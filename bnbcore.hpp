/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   bnbcore.hpp
 * Author: mposypkin
 *
 * Created on July 31, 2018, 4:46 PM
 */

#ifndef BNBCORE_HPP
#define BNBCORE_HPP

#include <algorithm>
#include <expression/expr.hpp>

#include "eqsystem.hpp"

using Box = std::vector<Interval<double>>;

std::ostream& operator<<(std::ostream& os, const Box& box) {
    for (auto& I : box) {
        std::cout << "[" << I.lb() << ":" << I.rb() << "]";
    }
    return os;
}

/**
 * A sort of a box
 */
enum class BoxSort {
    external,
    internal,
    boundary
};

std::ostream& operator<<(std::ostream& os, const BoxSort& obj) {
    switch (obj) {
        case BoxSort::boundary:
            os << "Boundary";
            break;
        case BoxSort::external:
            os << "External";
            break;
        case BoxSort::internal:
            os << "Internal";
            break;
        default:
            os << "Invalid";
    }
    return os;
}

/**
 * Computes the interval bounds
 * @param expr the evaluated expression 
 * @param box the box for paramters
 * @return resulting interval
 */
Interval<double> calcInterval(const Expr<Interval<double>> expr, const Box& box) {
    auto result = expr.calc(InterEvalAlg<double>(box));
    return result;
}

/**
 * Compute box max side (edge)
 * @param box the box to compute the max edge
 * @return 
 */
double maxEdge(const Box& box) {
    Interval<double> me = *std::max_element(box.begin(), box.end(),
            [](Interval<double> a, Interval<double> b) {
                return (a.rb() - a.lb()) < (b.rb() - b.lb());
            });

    return me.rb() - me.lb();
}

/**
 * Splits a box into tho subboxes and add them to the list of boxes
 * @param box the source box
 * @param boxlist the list of boxes
 */
void splitBox(const Box& box, std::vector<Box>& boxlist) {
    auto result = std::max_element(box.begin(), box.end(),
            [](Interval<double> a, Interval<double> b) {
                return (a.rb() - a.lb()) < (b.rb() - b.lb());
            });

    int i = result - box.begin();
    const double m = 0.5 * (box[i].lb() + box[i].rb());
    Box b1(box);
    b1[i] = Interval<double>(box[i].lb(), m);
    Box b2(box);
    b2[i] = Interval<double>(m, box[i].rb());
    boxlist.push_back(b1);
    boxlist.push_back(b2);
}

/**
 * Evaluates a sort of a box
 * @param es equality system
 * @param box a box to evaluate
 * @return a sort of a box
 */
BoxSort sortBox(const EqualitySystem& es, const Box& box) {
    for (auto ig : es.mIG) {
        Interval<double> I = calcInterval(ig, box);
        if ((I.lb() > 0) || (I.rb() < 0))
            return BoxSort::external;
    }
    return BoxSort::boundary;
}



#endif /* BNBCORE_HPP */

