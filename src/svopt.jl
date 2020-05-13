abstract type SVOptMethod end

struct SVHillClimb <: SVOptMethod end
struct SVGoldenSection <: SVOptMethod end

"Parameter used in SVGoldenSection method"
goldenRatio = (1 + sqrt(5)) / 2

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
Function optimizes single variable function using golden section method
(to find single minimum value)
"""
function (svhc::SVGoldenSection)(f, x0; ϵ = 1e-8, maxiter)
  # Get range based on provided starting point
  (a, b) = find_min_interval(f, x0)
  # Iterations
  i = 0
  while abs(a - b) >= ϵ && i ≤ maxiter
    i += 1
    # Calculate next 2 points
    x1 = b - (b - a) / goldenRatio
    x2 = a + (b - a) / goldenRatio
    # Compare function values from calculated points and determine where minimum is located
    if f(x1) < f(x2)
      b = x2
    else
      a = x1
    end
  end
  # Get the approximate value and return it
  result = (a + b) / 2
  resultY = f(result)
  return resultY, result
end
