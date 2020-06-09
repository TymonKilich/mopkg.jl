abstract type SVOptMethod end

struct mini <: SVOptMethod end

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

function (svhc::mini)(f, a, b; ϵ = 1e-5)
	ρ = (3-sqrt(5))/2
	w1 = a + (1.0 - ρ)*(b - a)
	w2 = a + ρ*(b - a)
	wL = min(w1, w2)
	wP = max(w1, w2)
	a1 = a
	b1 = b
	while abs(b1-a1) > ϵ
		if f(wL) < f(wP)
			b1 = wP
			wP = wL
			wL = b1 - ρ*(b1 - a1)
		else
			a1 = wL
			wL = wP
			wP = a1 + ρ*(b1 - a1)
		end
		wP1, wL1 = wP, wL
		wL = min(wP1, wL1)
		wP = max(wP1, wL1)
	end

return (a1+b1)/2.0, f((a1+b1)/2), a1, b1
		
end

