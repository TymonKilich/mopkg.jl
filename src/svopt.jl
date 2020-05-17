abstract type SVOptMethod end

struct SVHillClimb <: SVOptMethod end
struct IntervalHalving <: SVOptMethod end

"Finite central difference"
fdc(f, x; h=1e-5) = (f(x+h/2) - f(x-h/2))/h

"Second order central finite difference"
sfdc(f, x; h=1e-5) = (f(x+h) - 2f(x) + f(x-h))/h^2

"Find bracket with minimum"
function find_min_interval(f, x0; step=0.1, expandfactor=2.0)
    a, b = x0, x0 + step
    fa, fb = f(a), f(b)
    if fb > fa
        a, b = b, a
        fa, fb = fb, fa
        step = -step
    end
    while true
        c, fc = b + step, f(b + step)
        if fc > fb
            return a < c ? (a, c) : (c, a)
        end
        a, b = b, c
        fa, fb = fb, fc
        step = step*expandfactor
    end
end

function (svhc::SVHillClimb)(f, x0; ϵ, maxiter, dampingfactor=0.5, step=1.0)
    x, fp = x0, f(x0)
    s, fs = x0 + step, f(x0 + step)
    if fs > fp
        x, s = s, x
        fp, fs = fs, fp
        step = -step
    end
    i, fn = 0, Inf
    while abs(fn-fs) ≥ ϵ && i ≤ maxiter
        i += 1
        s += step
        fs = fn
        fn = f(s)
        if fn > fs
            step = -step*dampingfactor
        end
    end
    return fn, s
end

"""
Function optimizes single variable function using Interval Halving Method.
"""
function (svhc::IntervalHalving)(f, x0; ϵ, maxiter)
    a, b = find_min_interval(f, x0)
    x_m = (a + b) / 2
    l = b - a
    i = 0
    while abs(l) >= ϵ && i <= maxiter
        ym = f(x_m)
        x1 = a + (l/4)
        x2 = b - (l/4)
        y1 = f(x1)
        y2 = f(x2)
        if y1 < ym
            b = x_m
            x_m = x1
        elseif y2 < ym
            a = x_m
            x_m = x2
        else
            a = x1
            b = x2
        end
        l = b - a
        i += 1
    end
    return f(x_m), x_m
end


