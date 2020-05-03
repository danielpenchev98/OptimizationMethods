using LinearAlgebra
#β is momentum decay, v is momentum
function accelerated_steepest_descent(f,∇f,x;β=0.4,ϵ=0.005)
    x_prev, x_curr = NaN, x
    iteration = 1
    v = zeros(length(x))
    while iteration==1 || abs(f(x_curr)-f(x_prev))>ϵ
        direction = -∇f(x_curr)
        direction = direction/norm(direction)
        learning_rate = strong_backtracking_line_search(f,∇f,x_curr,direction)
        v = β*v + learning_rate * direction
        x_prev = x_curr
        x_curr = x_prev + v
        println(x_curr)
        iteration=iteration+1
    end
    println(iteration)
    return x_curr
end

using ReverseDiff
x0 = [-1.2,1]
rosenbrock(x; a=1.0, b=100.0) = (a-x[1])^2 + b*(x[2]-x[1]^2)^2
∇f = x -> ReverseDiff.gradient(rosenbrock,x)
println(accelerated_steepest_descent(rosenbrock,∇f,x0))
