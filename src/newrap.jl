struct SVNewtonRaphson <: SVOptMethod end

using Calculus

"Newton-Raphson method"
function (svhc::SVNewtonRaphson)(f, x; ϵ, maxiter)
	fx = f(x)
	i = 0
	d = Calculus.derivative(f, x)
	dd = Calculus.second_derivative(f, x)
	h = d/dd
	while abs(h) >= ϵ && i <= maxiter
		d = Calculus.derivative(f, x)
		dd = Calculus.second_derivative(f, x)
		h = d/dd
		x = x - h
		fx = f(x)
		i = i + 1
	end
	return fx, x
end

