abstract type SVOptMethod end

struct SVFiboSearch <: SVOptMethod end
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

function (svhc::SVFiboSearch)(f, a, b; eps)
    L = b - a
    k = 1
    #ilość potrzebnych wyliczeń wartości funkcji przy danej dokładności eps
    FN1 = trunc(Int, 2 * L / eps)
    i = 2
    F = [1, 1]
    while F[i] < FN1
        append!(F, F[i] + F[i-1])
        i = i + 1
    end
    N = size(F)[1] - 1
    while k < N
        l = F[N - k + 1] / F[N + 1]
        x1 = a + l
        x2 = b - l
        if f(x1) < f(x2)
            b = x2
        elseif f(x1) > f(x2)
            a = x1
        else
            a = x1
            b = x2
        end
        k = k + 1
    end
    
    return a, b
end
