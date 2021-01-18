#nonlinear conjugate gradient method

#Polak-Ribiere formula for β - it's possible to cycle indefinitely,
#in that case β will be reset to 0, deleting the data about the past search directions
#after about n itarations(where the dimension of the function is n) because we will
#run out of linear independent directions
function conjugate_gradients_method(f,∇f,x;ϵ=0.00005)
    path = []
    x_curr,x_prev, β, iteration = x,NaN, 0, 1
    push!(path,x_curr)
    direction = -∇f(x_curr)
    while iteration==1 || abs(f(x_curr) - f(x_prev)) > ϵ
        temp = f(x_curr)
        println("(x, y)=($x_curr, $temp)")
        direction = -∇f(x_curr) + β*direction
        α = strong_backtracking_line_search(f,∇f,x_curr,direction)
        x_prev = x_curr
        x_curr = x_prev + α*direction
        #The Wolfe conditions arent enough so β = max(0,formula)
        β = max(0, dot(∇f(x_curr),(∇f(x_curr)-∇f(x_prev)))/(∇f(x_prev)'*∇f(x_prev)))
        iteration=iteration+1
        push!(path,x_curr)
    end
    return path
end

using ReverseDiff
x0 = [7.0,20.0]
rosenbrock(x; a=1.0, b=100.0) = (a-x[1])^2 + b*(x[2]-x[1]^2)^2
branin(x; a=1, b=5.1/(4π^2), c=5/π, r=6, s=10, t=1/(8π)) = a*(x[2]-b*x[1]^2+c*x[1]-r)^2 + s*(1-t)*cos(x[1]) + s
∇f = x -> ReverseDiff.gradient(branin,x)


path = conjugate_gradients_method(branin,∇f,x0)
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
