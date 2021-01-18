using LinearAlgebra
function QuasiNewton(f,x,∇f, ϵ=0.00000001)
    path = []
    push!(path,x)
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
        push!(path,x_curr)
        δ ,γ = x_curr - x_prev, ∇f(x_curr) - ∇f(x_prev)

        if iteration == 1
            H = γ'*δ/(γ'*γ) * LinearAlgebra.I
        end
        H = BFGSFormula(H,δ,γ)

        direction = -H*∇f(x_curr)
        iteration = iteration + 1
    end
    println(iteration)
    return path
end

    DFPFormula(H,δ,γ) = H + (δ*δ')/(δ'*γ) - (H*γ*γ'*H)/(γ'*H*γ)

function BroydenFormula(H,δ,γ,ϕ)
    v = ((γ'*H*γ)^(1.0/2.0)) * (δ / (δ'*γ) - H*γ / (γ'*H*γ))
    return  DFPFormula(H,δ,γ) + ϕ*(v*v')
end

BFGSFormula(H,δ,γ) = BroydenFormula(H,δ,γ,1.0)

using ReverseDiff
x0 = [7.0,20.0]
rosenbrock(x; a=1.0, b=5.0) = (a-x[1])^2 + b*(x[2]-x[1]^2)^2
∇f = x -> ReverseDiff.gradient(rosenbrock,x)
branin(x; a=1, b=5.1/(4π^2), c=5/π, r=6, s=10, t=1/(8π)) = a*(x[2]-b*x[1]^2+c*x[1]-r)^2 + s*(1-t)*cos(x[1]) + s
∇f = x -> ReverseDiff.gradient(branin,x)


path = QuasiNewton(branin,x0,∇f)

using Plots
Plots.pyplot() # back to pyplot
x = -10:0.01:20.0
y = -10.0:0.01:20.0
g(x,y;a=1, b=5.1/(4π^2), c=5/π, r=6, s=10, t=1/(8π)) = begin
    a*(y-b*x^2+c*x-r)^2 + s*(1-t)*cos(x) + s
#        (1.0-x)^2 + 5.0*(y-x^2)^2
    end
X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(y, 1, length(x))
Z = map(g, X, Y)
default(size=(1000,600))

contour(x, y, g,levels=100,linewidth=2.0)
firstParam =[]
secondParam = []
push!(firstParam,path[1][1])
push!(secondParam,path[1][2])
for i in 2:size(path)[1]
    push!(firstParam,path[i][1])
    push!(secondParam,path[i][2])

end
scatter!(firstParam,secondParam, label="checkpoints")
plot!(firstParam,secondParam,label="path")
