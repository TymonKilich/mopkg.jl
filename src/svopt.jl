abstract type SVOptMethod end

struct IntervalHalving <: SVOptMethod end
struct SVHillClimb <: SVOptMethod end

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
Interval Halving Method for optimizing function.
"""

function (svhc::IntervalHalving)(f, x0; ϵ, maxiter)
    start, ending = find_min_interval(f, x0)
    middle = (start + ending) / 2
    i = 0
    diff = ending - start
    while abs(diff) >= ϵ && i <= maxiter
        x_1 = start  + (diff/4)
        x_2 = ending - (diff/4)
        y_1 = f(x_1)
        y_2 = f(x_2)
        y_middle = f(middle)
        if y_1 < y_middle
            ending = middle
            middle = x_1
        elseif y_2 < y_middle
            start = middle
            middle = x_2
        else
            start = x_1
            ending = x_2
        end
        diff = ending - start
        i += 1
    end
    return f(middle), middle
end





