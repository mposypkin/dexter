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

/**
 * Evaluates the <g(box)> which is a maximal deviation of constraint functions 
 * from zero
 * @param es equality system
 * @param box a box to evaluate
 * @return deviation
 */
double getAbs(const EqualitySystem& es, const Box& box) {
    double dev = 0;
    for (auto ig : es.mIG) {
        Interval<double> I = calcInterval(ig, box);
        dev = std::max(dev, std::abs(I.lb()));
        dev = std::max(dev, std::abs(I.rb()));
    }
    return dev;
}

/**
 * Computes the center of the box
 * @param box the box to find a center
 * @return the central point
 */
std::vector<double> getCenter(const Box& box) {
    std::vector<double> c;
    for (auto i : box) {
        const double ci = 0.5 * (i.lb() + i.rb());
        c.push_back(ci);
    }
    return c;
}

/**
 * Checks whether the box contains the 
 * @param es
 * @param box
 * @return true if zero was found, false otherwise
 */
bool containsRoot(const EqualitySystem& es, const Box& box) {
    constexpr double vzero = 1e-3;
    const double me = maxEdge(box);
    const double hmin = 1e-3 * me;
    auto convol = [ = ](const std::vector<double>& x){
        double v = 0;
        for (auto g : es.mG) {
            const double u = g.calc(FuncAlg<double>(x));
            v += u * u;
        }
        return v;
    };
    auto step = [&box, convol](std::vector<double>& x, double& v, int i, double h) {
        auto xn = x;
        xn[i] = x[i] + h;
        xn[i] = std::min(xn[i], box[i].rb());
        xn[i] = std::max(xn[i], box[i].lb());
        double vn = convol(xn);
        bool rv = false;
        if (vn < v) {
            v = vn;
            x[i] = xn[i];
            rv = true;
        }
        return rv;
    };

    auto dec = [](double&h) {
        h *= 0.5;
    };

    std::vector<double> x = getCenter(box);
    bool rv = true;
    double h = 0.25 * me;
    double v = convol(x);
    const int n = box.size();

    while ((h > hmin) && (v > vzero)) {
        for (int i = 0; i < n; i++) {
            bool rv = step(x, v, i, h) || step(x, v, i, -h);
            if (!rv) {
                dec(h);
            }
        }
    }
    if (v <= vzero)
        rv = true;
    else
        rv = false;
    //    std::cout << "# v = " << v;
    //    std::cout << ", rv = " << rv << "\n";
    return rv;
}

#endif /* BNBCORE_HPP */

