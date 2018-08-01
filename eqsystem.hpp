/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   eqsystem.hpp
 * Author: mposypkin
 *
 * Created on August 1, 2018, 4:46 PM
 */

#ifndef EQSYSTEM_HPP
#define EQSYSTEM_HPP

#include <vector>
#include <expression/expr.hpp>

using namespace snowgoose::expression;
using namespace snowgoose::interval;

struct EqualitySystem {
    std::vector<Expr<double>> mG;
    std::vector<Expr<Interval<double>>> mIG;
};

#endif /* EQSYSTEM_HPP */

