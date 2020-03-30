#TODO add restart optimization

function direct_search_simplex(f, x0; α=1.0, β=0.5, γ=2.0, ϵ=0.01, iter_max = 3200)
    simplex = generate_simplex(x0)
    iter = 1

    while iter < iter_max
        iter = iter + 1

        sort!(simplex, by=x -> f(x))
        x_worst = simplex[length(simplex)]
        x_second_worst = simplex[length(simplex)-1]
        x_best = simplex[1]

        centroid = (sum(simplex) - x_worst) / (length(simplex)-1)

        #if std(simplex)

        new_point = NaN

        #reflection phase
        x_reflection = centroid + α*(centroid-x_worst)

        #Expasion phase
        if f(x_reflection) < f(x_best)
            x_expasion = centroid + γ*(x_reflection - centroid)
            new_point = f(x_expasion)>f(x_reflection) ? x_expasion : x_reflection

        elseif f(x_reflection) <= f(x_second_worst)
            new_point = x_reflection

        #contraction
        else
            x_contaction = centroid + β*(x_worst - centroid)

            #Shrink contraction
            if f(x_contaction) > f(x_worst)
                simplex = shrink_contraction(simplex,x_best)
                continue
            end

            new_point = x_contaction
        end
        deleteat!(simplex,length(simplex))
        push!(simplex,new_point)
    end
    return simplex
end

function shrink_contraction(simplex, x_best; δ=0.5)
    shrinked_simplex = Array{Float64,1}[]
    for i in 1:length(x_best)+1
        push!(shrinked_simplex, x_best + δ*(simplex[i]-x_best))
    end
    return shrinked_simplex
end

function generate_simplex(x; α=0.05, β=0.00025)
    simplex = Array{Float64,1}[]
    push!(simplex,x)
    for i in 1:length(x)
        point = copy(x)
        point[i] = point[i] + ( x[i]==0 ? β : α )
        push!(simplex, point)
    end
    return simplex
end

x=[-1.2,1.0]
simplex = generate_simplex(x)
f(x)=sum(x)
println(sum(simplex))
println(sort(simplex, by=x -> f(x)))

rosenbrock(x; a=1.0, b=100) = (a-x[1])^2 + b*(x[2]-x[1]^2)^2
println(direct_search_simplex(rosenbrock,x))
