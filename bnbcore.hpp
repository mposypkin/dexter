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

#include <expression/expr.hpp>

Interval<double> calcInterval(const Expr<Interval<double>> expr, const std::vector<Interval<double>>& vars) {
    auto result = expr.calc(InterEvalAlg<double>(vars));
    return result;
}


#endif /* BNBCORE_HPP */

