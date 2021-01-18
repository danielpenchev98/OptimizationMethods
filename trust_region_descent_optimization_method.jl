function trust_region_descent(f,∇f,H,x,max_iterations; ϵ=0.000005, ƞ1=0.25, ƞ2=0.75, γ1=0.5,γ2=2.0, Δ=1.0, Δmax=10.0)
    path = []
    x_curr, y_curr, x_prev, y_prev = x, f(x), NaN, NaN
    push!(path,x_curr)
    iteration = 0
    while iteration <= max_iterations && ( iteration==0 || abs(y_curr-y_prev) > ϵ)
        x′, y′ = solve_trust_region_subproblem(f,∇f, H, x_curr, Δ)
        r = (y_curr - f(x′)) / (y_curr - y′)
        println("(x′, y′)=($x′, $y′)")
        if r < ƞ1
            Δ *= γ1
        else
            if r > ƞ2
                Δ = min(Δ*γ2,Δmax)
            end
            x_curr, y_curr, x_prev, y_prev = x′, y′, x_curr, y_curr

            δ, γ = x_curr-x_prev, ∇f(x_curr)-∇f(x_prev)

            H = updateHessian(H,δ,γ,iteration)
            push!(path,x_curr)
        end
        iteration+=1
    end
    return path
end

function updateHessian(H,δ,γ,iteration)
    if δ'*γ <= 0
        return H
    end

    if iteration == 1
        H = γ'*δ/(γ'*γ) * LinearAlgebra.I
    end
    H = BFSG(H,δ,γ)
end

BFSG(H,δ,γ) = H + (γ*γ')/(γ'*δ) - H*δ*transpose(H*δ)/(δ'*H*δ)

function solve_trust_region_subproblem(f,∇f, H, xk, Δ; ΔlargeThreshold = 8.0)
    m(p,xk) = f(xk) + transpose(p)*∇f(xk) + transpose(p)*H*p/2
    if Δ >= ΔlargeThreshold && transpose(∇f(xk))*∇f(xk)/(transpose(∇f(xk))*H*∇f(xk)) > 0
        println(Δ)
        pB = -inv(H)*∇f(xk)
        return xk + pB, m(pB,xk)
    else
        p = calculate_cauchy_point(∇f, H, xk, Δ)
        return xk + p, m(p,xk)
    end
end


function calculate_cauchy_point(∇f, H, xk, Δ)
    steepestDescentPoint = -Δ*∇f(xk)/norm(∇f(xk))
    step = 1
    positiveDefiniteTest = transpose(∇f(xk))*H*∇f(xk)
    if positiveDefiniteTest > 0
        step = min((norm(∇f(xk))^3) / (Δ * positiveDefiniteTest),1)
    end
    return step * steepestDescentPoint
end


using ReverseDiff
x0 = [-1.5,-1.5]
rosenbrock(x; a=1.0, b=5.0) = (a-x[1])^2 + b*(x[2]-x[1]^2)^2
∇f = x -> ReverseDiff.gradient(rosenbrock,x)
branin(x; a=1, b=5.1/(4π^2), c=5/π, r=6, s=10, t=1/(8π)) = a*(x[2]-b*x[1]^2+c*x[1]-r)^2 + s*(1-t)*cos(x[1]) + s
∇f = x -> ReverseDiff.gradient(branin,x)


path = trust_region_descent(rosenbrock, ∇f, LinearAlgebra.I , x0, 500)
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
