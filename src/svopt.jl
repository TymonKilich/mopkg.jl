abstract type SVOptMethod end

struct SVHillClimb <: SVOptMethod end
struct Goldss <: SVOptMethod end

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

function (svhc::Goldss)(f,x0;ϵ=1e-6,maxiter,τ=0.618)
	a=100-x0
	b=-a
    fx,fb=f(a),f(b)
    x1,x2=b-(τ)*(b-a),a+(τ)*(b-a)
    err=Inf
	i=0
    while err>ϵ && i ≤ maxiter
	    i +=1
        if f(x2)>f(x1)
            b=x2
            x2=x1
            x1=a+(1-τ)*(b-a)
        elseif f(x2)<f(x1)
            a=x1
            x1=x2
            x2=b-(1-τ)*(b-a)
        else
            a=(x1+x2)/2
            b=a
        end
        err=2*abs((b-a)/(b+a))
    end
    xsr=(x1+x2)/2
    return f(xsr),xsr
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
