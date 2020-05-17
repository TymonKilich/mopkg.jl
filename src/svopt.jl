abstract type SVOptMethod end

struct SVPowellMinimize <: SVOptMethod end

function p(f, x, v)

    a = v[1]
    b = v[2]
    c = v[3]

    fa = f(a)
    fb = f(b)
    fc = f(c)

    y = fa*(x-b)*(x-c)/((a-b)*(a-c))+fb*(x-a)*(x-c)/((b-a)*(b-c))+fc*(x-a)*(x-b)/((c-a)*(c-b))

    return y
end

function p_xm(f, v)

    a = v[1]
    b = v[2]
    c = v[3]

    fa = f(a)
    fb = f(b)
    fc = f(c)

    xm = 0.5*((b*b-c*c)*fa+(c*c-a*a)*fb+(a*a-b*b)*fc)/((b-c)*fa+(c-a)*fb+(a-b)*fc)

    return xm
end

function is_peak_min(f, xm, v)

    if(p(f, xm, v)<p(f, xm+0.1, v))
        return true
    else
        return false
    end
end

# PowellMinimize
# minimizing function
# arguments:
# f - function to be minimized
# a - start point
# h - initial step size
# ϵ - epsilon (accuracy)
# nmax - maximum steps
# returns value and localization of local minimum
function(svhc::SVPowellMinimize)(f, a; h = 0.5, ϵ = 1e-8, nmax = 100)

    b = a+h

    if(f(a)<f(b))
        c = a-h*2
    else
        c = a+3*h
    end

    v = [a, b, c]

    xm = p_xm(f, v)


     for i in 1:nmax

        is_minimum = is_peak_min(f, xm, v)

        println(i, " xm: ", xm, " is minimum ", is_minimum)

        for j in v
            if(abs(xm-j)<ϵ)
                return f(xm), xm
            end
        end

        if(is_minimum)
            k = 0
            old_dist = 0
            for j in 1:3
                dist = abs(xm-v[j])
                if(dist>old_dist)
                    k = j
                    old_dist = dist
                end
            end
            v[k] = xm
            xm_new = p_xm(f, v)
            if(p(f, xm_new, v) < p(f, xm, v))
                xm = xm_new
            else
                return f(xm), xm
            end
        else
            k = 0
            old_dist = 999
            for j in 1:3
                dist = abs(xm-v[j])
                if(dist<old_dist)
                    k = j
                    old_dist = dist
                end
            end
            if(xm>v[k])
                v[k] -= xm
            else
                v[k] += xm
            end

            xm_new = p_xm(f, v)
            xm = xm_new
        end
    end

    return f(xm), xm

end
