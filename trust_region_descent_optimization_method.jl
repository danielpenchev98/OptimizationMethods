function trust_region_descent(f,∇f,H,x,k_max; ƞ1=0.25, ƞ2=0.5, γ1=0.5,γ2=2.0, δ=1.0)
    y=f(x)
    for k in 1 : k_max
        x′, y′ = solve_trust_region_subproblem(∇f, H, x, δ)
        r = (y - f(x′)) / (y - y′)
        println("(x′, y′)=($x′, $y′)")
        if r < ƞ1
            δ *= γ1
        else
            x, y = x′, y′
            if r > ƞ2
                δ *= γ2
            end
        end
    end
    return x
end

using Convex, SCS
function solve_trust_region_subproblem(∇f, H, x0, δ)
    x = Variable(length(x0))
    p = Convex.minimize(dot(∇f(x0),(x-x0)) + Convex.quadform(x-x0, H(x0))/2)
    p.constraints += norm(x-x0) <= δ
    #there isnt a default solver - using SCS which is like a swiss knife
    Convex.solve!(p,SCS.Optimizer(verbose=false))
    return (x.value, p.optval)
end

using ReverseDiff

x0 = [-1.5,-1.5]
rosenbrock(x; a=1.0, b=5.0) = (a-x[1])^2 + b*(x[2]-x[1]^2)^2
∇f = x -> ReverseDiff.gradient(rosenbrock,x)
H = x -> ReverseDiff.hessian(rosenbrock,x)
println(trust_region_descent(rosenbrock, ∇f, H, x0, 10))
