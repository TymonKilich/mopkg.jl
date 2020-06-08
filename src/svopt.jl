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

function (svhc::IntervalHalving)(f, x0; ϵ, maxiter)
    s, e = find_min_interval(f, x0)
    mid = (s + e) / 2
    i = 0
    d = e - s
    while abs(d) >= ϵ && i <= maxiter
        x1 = s  + (d/4)
        x2 = e - (d/4)
        y1 = f(x1)
        y2 = f(x2)
        ym = f(mid)
        if y1 < ym
            e = mid
            mid = x1
        elseif y2 < ym
            s = mid
            mid = x2
        else
            s = x1
            e = x2
        end
        d = e - s
        i += 1
    end
    return f(mid), mid
end
