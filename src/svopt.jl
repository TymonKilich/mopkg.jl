abstract type SVOptMethod end

struct SV_Golden_Search <: SVOptMethod end

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


function (svhc::SV_Golden_Search)(f,x0;ϵ,maxiter,R=0.618,xL=50-x0,xU=-xL)
    w1 =  f(xL) 
    w2 =  f(xU)
    d = R * (xU-xL)
    t = 1 -R
    x1 = xU-d
    x2 = xL+d
	i=0
    while abs(w1 - w2)>ϵ && i ≤ maxiter
	    i +=1
        if f(x2)<f(x1)
            xL=x1
            x1=x2
            x2=xU- t*(xU-xL)
        elseif f(x2)>f(x1)
            xU=x2
            x2=x1
            x1=xL+ t*(xU-xL)
        else
            xU=xL
        end
        
    end
    return f(((x1+x2)/2)), (x1+x2)/2
end
