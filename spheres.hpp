/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   spheres.hpp
 * Author: mposypkin
 *
 * Created on August 2, 2018, 4:11 PM
 */

#ifndef SPHERES_HPP
#define SPHERES_HPP

#include <vector>
#include <expression/expr.hpp>
#include "eqsystem.hpp"

using namespace snowgoose::expression;
using namespace snowgoose::interval;

struct Spheres : public EqualitySystem {
    
    static constexpr double r1 = 5;
    static constexpr double r2 = 5;
    static constexpr double d = 8;
    

    Spheres() {
        EqualitySystem::mG.push_back(g1<double>());
        EqualitySystem::mG.push_back(g2<double>());
        EqualitySystem::mIG.push_back(g1<Interval<double>>());
        EqualitySystem::mIG.push_back(g2<Interval<double>>());
    }

    template <class T> Expr<T> g1() {
        Expr<T> x;
        Expr<T> y = sqr(x[0]) + sqr(x[1]) + sqr(x[2]) - r1 * r1;
        return y;
    }

    template <class T> Expr<T> g2() {
        Expr<T> x;
        Expr<T> y = sqr(x[0]) + sqr(x[1]) + sqr(x[2] - d) - r2 * r2;
        return y;
    }


};




#endif /* SPHERES_HPP */

