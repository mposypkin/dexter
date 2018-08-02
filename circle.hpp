/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   circle.hpp
 * Author: mposypkin
 *
 * Created on August 2, 2018, 3:16 PM
 */

#ifndef CIRCLE_HPP
#define CIRCLE_HPP

#include <vector>
#include <expression/expr.hpp>
#include "eqsystem.hpp"

using namespace snowgoose::expression;
using namespace snowgoose::interval;

struct Circle : public EqualitySystem {
    // xp = x[0], yp = x[1], q1 = x[2], q2 = x[3], q3 = x[4], q4 = x[5]

    static constexpr double r = 3;


    Circle() {
        EqualitySystem::mG.push_back(g1<double>());
        EqualitySystem::mIG.push_back(g1<Interval<double>>());
    }
    
    template <class T> Expr<T> g1() {
        Expr<T> x;
        Expr<T> y = sqr(x[0]) + sqr(x[1]) - r * r;
        return y;
    }

};


#endif /* CIRCLE_HPP */

