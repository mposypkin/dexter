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
#include <src/optimization.h>

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

static const EqualitySystem* pes;

/**
 * Checks whether the box contains the roots by coordinate descent
 * @param es equation system
 * @param box box constraints
 * @return 
 */
bool containsRootNewt(const EqualitySystem& es, const Box& box) {
    std::vector<double> c = getCenter(box);
    pes = &es;
    const int n = box.size();
    const int m = es.mG.size();
    //    for (auto dg : es.mDG) {
    //        std::cout << "OK in\n";
    //        ValDer<double> der = dg.calc(ValDerAlg<double>(x));
    //        Grad<double> grad = der.grad();
    //        for (int i = 0; i < n; i++) {
    //            std::cout << "a [" << i << "] = " << grad.getGrad()[i] << "\n";
    //        }
    //
    //    }

    auto function1_fvec = [](const alglib::real_1d_array &x, alglib::real_1d_array &fi, void *ptr) {
        //
        // this callback calculates
        // f0(x0,x1) = 100*(x0+3)^4,
        // f1(x0,x1) = (x1-3)^4
        //
        //        fi[0] = x[0] * x[0] + x[1] * x[1] + (x[2] - 2) * (x[2] - 2) - 9;
        //        fi[1] = x[0] * x[0] + x[1] * x[1] + (x[2] + 2) * (x[2] + 2) - 9;
        //fi[0] = fi[0] * fi[0]; 
        //fi[0] = fi[1] * fi[1]; 
        std::vector<double> xx(x.getcontent(), x.getcontent() + x.length());
        int i = 0;
        for (auto g : pes->mG) {
            const double v = g.calc(FuncAlg<double>(xx));
            fi[i++] = v;
        }
    };

    auto function1_jac = [](const alglib::real_1d_array &x, alglib::real_1d_array &fi, alglib::real_2d_array &jac, void *ptr) {
        //
        // this callback calculates
        // f0(x0,x1) = 100*(x0+3)^4,
        // f1(x0,x1) = (x1-3)^4
        // and Jacobian matrix J = [dfi/dxj]
        //
        std::vector<double> xx(x.getcontent(), x.getcontent() + x.length());
        int i = 0;
        const int n = x.length();
        //        std::cout << "n = " << xx.size() << "\n";
        //        std::cout.flush();
        for (auto dg : pes->mDG) {
            ValDer<double> der = dg.calc(ValDerAlg<double>(xx));
            const double v = der.value();
            fi[i] = v;
            Grad<double> grad = der.grad();
            for (int j = 0; j < n; j++) {
                jac[i][j] = grad.getGrad()[j];
            }
            i++;
        }

        //        jac[0][0] = 2 * x[0];
        //        jac[0][1] = 2 * x[1];
        //        jac[0][2] = 2 * x[2];
        //        jac[1][0] = 2 * x[0];
        //        jac[1][1] = 2 * x[1];
        //        jac[1][2] = 2 * x[2];
    };

    alglib::real_1d_array x;
    x.setlength(n);
    alglib::real_1d_array bndl;
    bndl.setlength(n);
    alglib::real_1d_array bndu;
    bndu.setlength(n);
    for (int i = 0; i < n; i++) {
        x[i] = c[i];
        bndl[i] = box[i].lb();
        bndu[i] = box[i].rb();
    }
    double epsx = 0.0000000001;

    alglib::minlmstate state;
    alglib::ae_int_t maxits = 0;
    alglib::minlmreport rep;
    //alglib::minlmcreatev(m, x, 0.0001, state);
    alglib::minlmcreatevj(m, x, state);
    minlmsetcond(state, epsx, maxits);
    minlmsetbc(state, bndl, bndu);

    //alglib::minlmoptimize(state, function1_fvec);
    alglib::minlmoptimize(state, function1_fvec, function1_jac);
    alglib::minlmresults(state, x, rep);
    //    std::cout << x.tostring(8) << "\n";
    std::vector<double> xx(x.getcontent(), x.getcontent() + x.length());
    double s = 0;
    for (auto g : pes->mG) {
        //std::cout << g.calc(FuncAlg<double>(xx)) << "\n";
        double ts = g.calc(FuncAlg<double>(xx));
        s += ts * ts;
    }
    bool good = true;
//    std::cout << "s = " << s << "\n";
    if (s > 1e-4) {
        good = false;
    } else {
        //        for (int i = 0; i < n; i++) {
        //            if ((xx[i] < box[i].lb()) || (xx[i] > box[i].rb())) {
        //                good = false;
        //                break;
        //            }
        //        }
    }
    return good;
}

/**
 * Checks whether the box contains the roots by coordinate descent
 * @param es equation system
 * @param box box to search
 * @return true if zero was found, false otherwise
 */
bool containsRootCdesc(const EqualitySystem& es, const Box& box) {
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

/**
 * Checks whether a box contains a root of a system
 * @param es equations system
 * @param box box to search
 * @return 
 */
bool containsRoot(const EqualitySystem& es, const Box& box) {
    return containsRootNewt(es, box);
//        return containsRootCdesc(es, box);
}

#endif /* BNBCORE_HPP */

