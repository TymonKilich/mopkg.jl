struct SVPowell <: SVOptMethod end

"Function returns minumum of a function with Powell's interpolation method"

function (svpl::SVPowell)(f, x0; ϵ=1e-5, maxiter = 10^6, step=1.0)

    c = x0 + step

    if f(x0)<f(c)
        a = x0 - step
    else
        a = x0 + 2*step
        a, c = c, a
    end

    b = x0

    maxStep = 5.0

    i = 1

    while i < maxiter

        first = (1/2)
        second = (b^2-c^2) * f(a) + (c^2-a^2) * f(b) + (a^2-b^2)f(c)
        third = (b-c) * f(a) + (c-a) * f(b) + (a-b) * f(c)

        d = first*second/third

        first = (b-c) * f(a) + (c-a) * f(b) + (a-b) * f(c)
        second = (a-b) * (b-c) * (c-a)

        if first/second < 0
        
            if abs(a-d)<ϵ 
                return (f(a), a)
            end

            if abs(b-d)<ϵ 
                return (f(b), b)
            end

            if abs(c-d)<ϵ 
                return (f(c), c)
            end

            list = sort([(a,f(a)), (b,f(b)), (c,f(c))], by = x -> x[2])


            if f(d) > f(a)
                far = sort([(a,abs(a-d)), (b,abs(b-d)), (c,abs(c-d))], by = x -> x[2])
                final = sort([far[1][1], far[2][1], d])
                a = final[1]
                b = final[2]
                c = final[3]
            else
                a = d
                b = list[1][1]
                c = list[2][1]
            end
        else
            x0 = x0 + maxStep

            c = x0 + step

            if f(x0)<f(c)
                a = x0 - step
            else
                a = x0 + 2*step
                a, c = c, a
            end

            b = x0
        end
        i+=1
    end

    throw("Couldn't find minimum")

end

