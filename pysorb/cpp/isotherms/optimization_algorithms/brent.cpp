
#include "brent.h"
#include <iostream>

double MAX_ITER = 1500;

/**
 * Returns the first value in the same unit as de second one
 *
 * @param value1 Value to change sign.
 * @param value2 Value to find sign.
 * @return value1 in the same unit as value2
 */
double
sign(double value1, double value2)
{
    if (value2 >= 0)
        return fabs(value1);
    else
        return -fabs(value1);
}

std::tuple<double, double, double, double, double, double> bracket(std::function<double(double)> fun)
{
    double denom, fw;
    double grow_limit = 110.0;
    int maxiter = 1000;
    double xa = 0.0;
    double xb = 1.0;
    double _gold = 1.618034; // golden ratio: (1.0+sqrt(5.0))/2.0
    double _verysmall_num = 1e-21;
    double fa = fun(xa);
    double fb = fun(xb);
    double tmpx;
    double tmpf;
    if (fa < fb)
    {
        tmpx = xa;
        xa = xb;
        xb = tmpx;
        tmpf = fa;
        fa = fb;
        fb = tmpf;
    }
    double xc = xb + _gold * (xb - xa);
    double fc = fun(xc);
    int funcalls = 3;
    int iter = 0;
    while (fc < fb)
    {
        double tmp1 = (xb - xa) * (fb - fc);
        double tmp2 = (xb - xc) * (fb - fa);
        double val = tmp2 - tmp1;
        if (std::fabs(val) < _verysmall_num)
        {
            denom = 2.0 * _verysmall_num;
        }
        else
        {
            denom = 2.0 * val;
        }
        double w = xb - ((xb - xc) * tmp2 - (xb - xa) * tmp1) / denom;
        double wlim = xb + grow_limit * (xc - xb);
        if (iter > maxiter)
        {
            break;
        }
        iter += 1;
        if ((w - xc) * (xb - w) > 0.0)
        {
            fw = fun(w);
            funcalls++;
            if (fw < fc)
            {
                xa = xb;
                xb = w;
                fa = fb;
                fb = fw;

                return std::make_tuple(xa, xb, xc, fa, fb, fc);
            }
            else if (fw > fb)
            {
                xc = w;
                fc = fw;
                return std::make_tuple(xa, xb, xc, fa, fb, fc);
            }
            w = xc + _gold * (xc - xb);
            fw = fun(w);
            funcalls++;
        }
        else if ((w - wlim) * (wlim - xc) >= 0.0)
        {
            w = wlim;
            fw = fun(w);
            funcalls++;
        }
        else if ((w - wlim) * (xc - w) > 0.0)
        {
            fw = fun(w);
            funcalls++;

            if (fw < fc)
            {
                xb = xc;
                xc = w;
                w = xc + _gold * (xc - xb);
                fb = fc;
                fc = fw;
                fw = fun(w);
                funcalls++;
            }
        }
        else
        {
            w = xc + _gold * (xc - xb);
            fw = fun(w);
            funcalls++;
        }
        xa = xb;
        xb = xc;
        xc = w;
        fa = fb;
        fb = fc;
        fc = fw;
    }
    return std::make_tuple(xa, xb, xc, fa, fb, fc);
}

std::tuple<double, double, double, double> get_bracket_info(std::function<double(double)> fun, double x)
{
    double a, b, c, d, e, fa, fb, fc, fx, dx;

    fx = fun(x);

    if (x != 0.0)
        dx = x / 100.0;
    else
        dx = 1.0 / 100.0;

    // assign initial values
    a = x;
    fa = fx;
    b = x;
    fb = fx;

    // Finds a sign change

    int max_iter = 1000;
    int iter = 0;
    while ((fa > 0.) == (fb > 0.))
    {
        dx = 1.414213562373095048801 * dx;
        a = x - dx;
        fa = fun(a);

        if ((fa > 0) != (fb > 0) || iter == max_iter)
            break;

        b = x + dx;
        fb = fun(b);

        iter++;
    }

    // Initialization

    fa = fun(a);
    fb = fun(b);

    return std::make_tuple(a, b, fa, fb);
}

double brent_zeroin2(std::function<double(double)> fun, double x, double tol)
{
    double a, b, c, d, e, fa, fb, fc, tol1, xm, p, q, r, s, eps, fx, dx, tol2, xmid, rat, xa, xb, xc;
    double tmp1, tmp2, fv, fw, w, v, dx_temp, u, fu;
    eps = 1.1e-11;
    // std::tie(xa, xb, xc, fa, fb, fc) = bracket(fun);

    std::tie(a, b, fa, fb) = get_bracket_info(fun, x);
    fx = fun(x);
    w = v = x;
    fw = fv = fx;

    int iter = 0;
    int funcalls = 0;

    double deltax = 0.0;
    double _cg = 0.3819660;

    while (iter < MAX_ITER)
    {

        tol1 = tol * std::fabs(x) + eps;
        tol2 = 2.0 * tol1;

        xmid = 0.5 * (a + b);

        if (std::fabs(x - xmid) < (tol - 0.5 * (b - a)))
            break;

        if (std::fabs(deltax) <= tol1)
        {
            if (x >= xmid)
                deltax = a - x;
            else
                deltax = b - x;
            rat = _cg * deltax;
        }
        else
        {
            tmp1 = (x - w) * (fx - fv);
            tmp2 = (x - v) * (fx - fw);
            p = (x - v) * tmp2 - (x - w) * tmp1;
            tmp2 = 2.0 * (tmp2 - tmp1);
            if (tmp2 > 0)
                p = -p;
            tmp2 = std::fabs(tmp2);
            dx_temp = deltax;
            deltax = rat;

            if (((p > tmp2 * (a - x)) && (p < tmp2 * (b - x))) && (std::fabs(p) < std::fabs(0.5 * tmp2 * dx_temp)))
            {
                rat = p * 1.0 / tmp2;
                u = x + rat;
                if ((u - a) < tol2 || (b - u) < tol2)
                {
                    if (xmid - x >= 0)
                        rat = tol1;
                    else
                        rat = -tol1;
                }
                else
                {
                    if (x >= xmid)
                        deltax = a - x;
                    else
                        deltax = b - x;
                    rat = _cg * deltax;
                }
            }
        }
        if (std::fabs(rat) < tol1)
        {
            if (rat >= 0)
                u = x + tol1;
            else
                u = x - tol1;
        }
        else
        {
            u = x + rat;
        }

        fu = fun(u);
        funcalls++;

        if (fu > fx)
        {
            if (u < x)
                a = u;
            else
                b = u;

            if ((fu <= fw) || (w == x))
            {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            }
            else if ((fu <= fv) || (v == x) || (v == w))
            {
                v = u;
                fv = fu;
            }
        }
        else
        {
            if (u >= x)
                a = x;
            else
                b = x;

            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
        }
        iter++;
    }

    return x;
}

/**
 * Perform the Richard Brent's improvements to Dekker's zeroin algorithm
 *
 * @param fun Pointer to Function that receives a double and returns a double (function to be minimized).
 * @param x Initial estimate.
 * @param tol Minimum `TOL` .
 * @return Value of the value with minimal :fun: value.
 */
double brent_zeroin(std::function<double(double)> fun, double x, double tol)
{
    double a, b, c, d, e, fa, fb, fc, tol1, xm, p, q, r, s, fx, dx;

    std::tie(a, b, fa, fb) = get_bracket_info(fun, x);
    double eps = 1e-11;
    int iter = 0;

// Begin step
begin_20:
    c = a;
    fc = fa;
    d = b - a;
    e = d;
begin_30:
    if (std::fabs(fc) >= std::fabs(fb))
        goto convergence_test;
    a = b;
    b = c;
    c = a;
    fa = fb;
    fb = fc;
    fc = fa;

convergence_test:
    tol1 = 2.0 * eps * std::fabs(b) + 0.5 * tol;
    xm = 0.5 * (c - b);

    if (std::fabs(xm) <= tol1 || fb == 0.0)
        goto done;
    if (iter == MAX_ITER)
    {
        // std::cout << "Reached maximum iteractions: " << MAX_ITER << std::endl;
        goto done;
    }

    iter++;
    // Is bisection necessary
    if (std::fabs(e) < tol1)
        goto bisection;
    if (std::fabs(fa) <= std::fabs(fb))
        goto bisection;

    // Is quadratic size_terpolation possible
    if (a != c)
        goto inverse_quadratic_size_terpolation;

    s = fb / fa;
    p = 2.0 * xm * s;
    q = 1.0 - s;

    goto adjust_signs;

inverse_quadratic_size_terpolation:
    q = fa / fc;
    r = fb / fc;
    s = fb / fa;
    p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
    q = (q - 1.0) * (r - 1.0) * (s - 1.0);

adjust_signs:
    if (p > 0.0)
    {
        q = -q;
    }
    p = std::fabs(p);

    // Is size_terpolation acceptable
    if (2.0 * p >= (3.0 * xm * q - fabs(tol1 * q)) || p >= fabs(0.5 * e * q))
        goto bisection;
    e = d;
    d = p / q;
    goto complete_step;

bisection:
    d = xm;
    e = d;

complete_step:
    a = b;
    fa = fb;
    if (std::fabs(d) > tol1)
        b = b + d;
    else
        b = b + sign(tol1, xm);
    fb = fun(b);
    if ((fb * (fc / fabs(fc))) > 0.0)
        goto begin_20;
    goto begin_30;

done:
    return b;
}