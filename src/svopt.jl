abstract type SVOptMethod end

struct SVHillClimb <: SVOptMethod end
struct SVDivInHalf <: SVOptMethod end

"Finite central difference"
fdc(f, x; h=1e-5) = (f(x+h/2) - f(x-h/2))/h

"Second order central finite difference"
sfdc(f, x; h=1e-5) = (f(x+h) - 2f(x) + f(x-h))/h^2

"Find bracket with minimum"
function find_min_interval(f, x0; step=1, expandfactor=2.0)
#changed step value to integer as approximated interval
#is much better with integer boundries
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

function (svhc::SVDivInHalf)(f, x0; ϵ, maxiter)

	err = ϵ
	a, b = find_min_interval(f, x0)
	middle = (a+b)/2
	diff = b - a
	iter= 0

	while(abs(diff) >= err && iter < maxiter)
		iter += 1
		left = a + diff/4
		right = b - diff/4
		if(f(left) < f(middle))
			b = middle
			middle = left
		elseif(f(right) < f(middle))
			a = middle
			middle = right
		else
			a = left
			b = right
		end
		diff = b - a
	end
	return (f(middle), middle)
end
