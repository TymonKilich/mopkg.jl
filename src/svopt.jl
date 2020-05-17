abstract type SVOptMethod end

struct SVHillClimb <: SVOptMethod end
struct SVFiboMethod <: SVOptMethod end
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
function fibonacci_num(n) #generate fibonacci n-th number (recursive)
    return n < 2 ? n : fibonacci_num(n-1) + fibonacci_num(n-2)
end # function fibonacci_num

"""
    Method that find minimum of one argument function
    Arguments:
        f - function that we like to find minimum (ex. f(x) = x^2 - 2)
        a,b - interval beetwen we are looking for result
        tol - result tolerance
"""
function (svfb::SVFiboMethod)(f,a,b,tol=0.1)
    # 1 step - fiding n (number of iterations)
    n = 2
    nval = nval = (b-a)/fibonacci_num(n) # helper value for n calculation
    xres = 0 # result of calculation (x point)
    while nval >= (2*tol)
        nval = (b-a)/fibonacci_num(n)
        n+=1
    end # while
    while true
        # caclulate x1,x2 initial values
        x1 = b - ((fibonacci_num(n-1) * (b-a))/fibonacci_num(n))
        x2 = a + ((fibonacci_num(n-1) * (b-a))/fibonacci_num(n))
        # calculate values of function in x1,x2
        f1 = f(x1)
        f2 = f(x2)
        # eliminating areas
        if f1 < f2
            b = x2
            x2 = x1
            n = n - 1
            x1 = b - ((fibonacci_num(n-1) * (b-a))/fibonacci_num(n))
            xres = x1
        else
            a = x1
            x1 = x2
            n = n - 1
            x2 = a + ((fibonacci_num(n-1) * (b-a))/fibonacci_num(n))
            xres = x2
        end #if
        if abs(x2-x1) <= tol || n <= 2
            yres = f(xres)
            result = (xres,yres) # result as tuple
            return result # end of our function
        end # if
    end # while

end # function fibonacci
