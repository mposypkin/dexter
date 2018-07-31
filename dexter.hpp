/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   dexter.hpp
 * Author: mposypkin
 *
 * Created on July 31, 2018, 3:44 PM
 */

#ifndef DEXTER_HPP
#define DEXTER_HPP

#include <vector>
#include <expression/expr.hpp>

using namespace snowgoose::expression;
using namespace snowgoose::interval;

struct DexterKinematicEquations {
    // xp = x[0], yp = x[1], q1 = x[2], q2 = x[3], q3 = x[4], q4 = x[5]

    static constexpr double la = 1;

    static constexpr double lb = 1;

    static constexpr double d2 = 0.5;

    DexterKinematicEquations() {
        mG.push_back(g1<double>());
        mG.push_back(g2<double>());
        mG.push_back(g3<double>());
        mG.push_back(g4<double>());
        mIG.push_back(g1<Interval<double>>());
        mIG.push_back(g2<Interval<double>>());
        mIG.push_back(g3<Interval<double>>());
        mIG.push_back(g4<Interval<double>>());
    }
    
    template <class T> Expr<T> g1() {
        Expr<T> x;
        Expr<T> y = x[0] - la * cos(x[2]) - lb * cos(x[3]) - d2;
        return y;
    }

    template <class T> Expr<T> g2() {
        Expr<T> x;
        Expr<T> y = x[0] - la * cos(x[4]) - lb * cos(x[5]) + d2;
        return y;
    }

    template <class T> Expr<T> g3() {
        Expr<T> x;
        Expr<T> y = x[1] - la * sin(x[2]) - lb * sin(x[3]);
        return y;
    }

    template <class T> Expr<T> g4() {
        Expr<T> x;
        Expr<T> y = x[1] - la * sin(x[4]) - lb * sin(x[5]);
        return y;
    }

    std::vector<Expr<double>> mG;
    std::vector<Expr<Interval<double>>> mIG;
};

#endif /* DEXTER_HPP */

