abstract type SVOptMethod end

struct SVInterpol_DSC <: SVOptMethod end

function (inpol_dsc_::SVInterpol_DSC)(f, x0k; e, s)

    function rec_x(f, xp1k, s, p, n=1, xp2k=nothing)
        xnk = xp1k + (2 ^ (n - 1) * p * s)
        if f(xnk) > f(xp1k)
            return n, xp1k, xnk, xp2k
        end

        return rec_x(f, xnk, s, p, n + 1, xp1k)
    end

    function inpol_dsc(f, x0k, e, s, K=0.1)
        xp1k = x0k - s
        xn1k = x0k + s
        f0k = f(x0k)
        fn1k = f(xn1k)
        if f0k > fn1k
            p = 1
            n, xp1k, xnk, xp2k = rec_x(f, x0k, s, p)

        else
            fp1k = f(xp1k)
            if fp1k < f0k
                p = -1
                n, xp1k, xnk, xp2k = rec_x(f, x0k, s, p)
            
            elseif fp1k >= f0k && f0k <= fn1k
                x0kn1 = x0k + ((s * (fp1k - fn1k)) / (2 * (fp1k - (2 * f0k) + fn1k)))

                if s <= e
                    return x0kn1, f(x0kn1)
                else
                    inpol_dsc(f, x0kn1, e ,s * K)
                end
            end
        end

        xmk = xp1k + (2 ^ (n - 2) * p * s)
        fmk = f(xmk)

        if fmk >= f(xp1k)
            x0kn1 = xp1k + ((2 ^ (n - 2) * p * s * (f(xp2k) - fmk)) / (2 * (
                    f(xp2k) - 2 * f(xp1k) + fmk)))
        elseif fmk < f(xp1k)
            x0kn1 = xmk + (((2 ^ (n - 2)) * p * s * (f(xp1k) - f(xnk))) / (2 * (f(xp1k) - (2 * fmk) + f(xnk))))
        end

        if (2 ^ (n - 2)) * s <= e
            return x0kn1, f(x0kn1)
        else
            inpol_dsc(f, x0kn1, e,s * K)
        end
    end

    inpol_dsc(f, x0k, e, s)

end