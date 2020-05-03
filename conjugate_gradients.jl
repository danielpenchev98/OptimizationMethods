#nonlinear conjugate gradient method

#Polak-Ribiere formula for β - it's possible to cycle indefinitely,
#in that case β will be reset to 0, deleting the data about the past search directions
#after about n itarations(where the dimension of the function is n) because we will
#run out of linear independent directions
function conjugate_gradients_method(f,∇f,x;ϵ=0.005)
    x_curr,x_prev, β, iteration = x,NaN, 0, 1
    direction = -∇f(x_curr)
    while iteration==1 || abs(f(x_curr) - f(x_prev)) > ϵ
        direction = -∇f(x_curr) + β*direction
        alpha = strong_backtracking_line_search(f,∇f,x_curr,direction)
        x_prev = x_curr
        x_curr = x_prev + alpha*direction
        #The Wolfe conditions arent enough so β = max(0,formula)
        β = max(0, dot(∇f(x_curr),(∇f(x_curr)-∇f(x_prev)))/(∇f(x_prev)'*∇f(x_prev)))
        iteration=iteration+1
            println(iteration)
    end
    println(iteration)
    return x_curr
end

using ReverseDiff
x0 = [-1.2,1]
rosenbrock(x; a=1.0, b=100.0) = (a-x[1])^2 + b*(x[2]-x[1]^2)^2
∇f = x -> ReverseDiff.gradient(rosenbrock,x)
println(conjugate_gradients_method(rosenbrock,∇f,x0))
