using LinearAlgebra
function QuasiNewton(f,x,∇f, ϵ=0.00000001)
    #H - inverse Hessian approximation, γ = ∇f(x(i+1)) - ∇f(x(i)), δ = x(i+1) - x(i)
    H, δ, γ = LinearAlgebra.I, NaN, NaN;
    direction = -H*∇f(x)
    x_prev, x_curr, iteration = NaN, x, 1
    while iteration==1 || norm(∇f(x_curr)) > ϵ

        # strong Wolfe conditions guarantee the correctness of curvature condition δ'*γ>0
        # if curvature condition is fulfilled then Hessian is positive definite matrix
        alpha = strong_backtracking_line_search(f,∇f,x_curr,direction)
        x_prev = x_curr
        x_curr = x_prev + alpha*direction
        δ ,γ = x_curr - x_prev, ∇f(x_curr) - ∇f(x_prev)

        if iteration == 1
            H = γ'*δ/(γ'*γ) * LinearAlgebra.I
        end
        H = BFGSFormula(H,δ,γ)

        direction = -H*∇f(x_curr)
        iteration = iteration + 1
        println(iteration)
    end
    println(x_curr)
end

DFPFormula(H,δ,γ) = H + (δ*δ')/(δ'*γ) - (H*γ*γ'*H)/(γ'*H*γ)

function BroydenFormula(H,δ,γ,ϕ)
    v = ((γ'*H*γ)^(1.0/2.0)) * (δ / (δ'*γ) - H*γ / (γ'*H*γ))
    return  DFPFormula(H,δ,γ) + ϕ*(v*v')
end

BFGSFormula(H,δ,γ) = BroydenFormula(H,δ,γ,1.0)

using ReverseDiff
x0 = [-1.2,1]
rosenbrock(x; a=1.0, b=100.0) = (a-x[1])^2 + b*(x[2]-x[1]^2)^2
∇f = x -> ReverseDiff.gradient(rosenbrock,x)
QuasiNewton(rosenbrock,x0,∇f)
