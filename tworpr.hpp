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

#ifndef TWORPR_HPP
#define TWORPR_HPP

#include <vector>
#include <expression/expr.hpp>
#include "eqsystem.hpp"

using namespace snowgoose::expression;
using namespace snowgoose::interval;

struct TwoRPRKinematicEquations : public EqualitySystem {

    static constexpr double d2 = 3;

    TwoRPRKinematicEquations() {
        EqualitySystem::mG.push_back(g1<double>());
        EqualitySystem::mG.push_back(g2<double>());
        EqualitySystem::mIG.push_back(g1<Interval<double>>());
        EqualitySystem::mIG.push_back(g2<Interval<double>>());
    }
    
    template <class T> Expr<T> g1() {
        Expr<T> x;
        Expr<T> y = sqr(x[0] - d2) + sqr(x[1]) - sqr(x[2]);
        return y;
    }

    template <class T> Expr<T> g2() {
        Expr<T> x;
        Expr<T> y = sqr(x[0] + d2) + sqr(x[1]) - sqr(x[3]);
        return y;
    }

};

#endif /* DEXTER_HPP */

