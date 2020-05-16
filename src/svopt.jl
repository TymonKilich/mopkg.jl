abstract type SVOptMethod end

struct SVHillClimb <: SVOptMethod end
struct SVGoldenRatio <: SVOptMethod end

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

"function optimize single variable function with golden section"

function (goldenratio::SVGoldenRatio)(f, x0;
                  ϵ=1e-5, maxiter=100)
  golrat = 1.61803398875
  i = 1
  a, b = find_min_interval(f, x0)
  fa = f(a)
  fb = f(b)
  temp = golrat * (a-b)
  rev = 1 - golrat
  while abs(b - a) > ϵ && i ≤ maxiter
	#x1 =b-(golrat)*(b-a) 
	#x2 = a+(golrat)*(b-a)
    x1 = b-(b-a)/ golrat # + temp
	x2 = a+(b-a)/ golrat # - temp
	fx1 = f(x1)
	fx2 = f(x2)
    if fx1 > fx2
      a = x1
	  x1 = b-(b - a)/ golrat
	  fx1 = f(x1)
	  fa = f(a)
    elseif fx1 < fx2
      b = x2
	  x2 = a+(b - a)/ golrat
	  fx2 = f(x2)
	  fb = f(b)
    else
	a = (x1+x2)/2
	fa = f(a)
	b = a
	fb = f(b)
	end
	i  = i + 1
  end
  return f(((a+b)/2)), ((a+b)/2)
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

