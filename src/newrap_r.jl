struct SVNewtonRaphson_r <: SVOptMethod end
using Calculus

"Newton Raphson method"


function (svnr_r::SVNewtonRaphson_r)(f, x0; ϵ=1e-5, maxiter=1e6)
           x_new, x_old = x0, Inf
           f_new, f_old = Calculus.derivative(f, x_new), Inf

           i = 0

           while i <= maxiter && abs(x_new-x_old) >= ϵ && abs(f_new-f_old) >= ϵ
               x_pocket = x_new - Calculus.derivative(f, x_new)/Calculus.second_derivative(f,  x_new)
               x_new, x_old = x_pocket, x_new
               f_new, f_old = Calculus.derivative(f, x_new), f_new
               i += 1

           end

           if i >= maxiter
               error("Function itered over ",maxiter," times.")
           else
               return f(x_new) ,x_new
           end
end

