using LinearAlgebra
function direct_search_jeeve_hooke(f, x ,Δ0 ,τ ; α=2.0, ϵ=0.01)
    x_begin = x
    Δ = copy(Δ0)
    while true

        #Exploration seacth step
        x_curr,success = exploration(f,x_begin,Δ)
        if !success
            Δ = Δ/2
            for i in 1:length(Δ)
                if Δ[i] < τ[i]
                    return x_curr
                end
            end
            continue
        else
            Δ = Δ0
        end

       # Pattern move step
       while true
            x_tentative = α*x_curr - x_begin

            x_final,success = exploration(f,x_tentative,Δ)
            x_begin = x_curr

            #should we continue to search in that direction
            if f(x_final) >= f(x_curr)
                break
            end
            x_curr = x_final
        end

    end

    return x_begin
end

#search for better point
function exploration(f,x0,Δ)
    x_curr = x0
    success = false
    for i in 1:length(x_curr)

        x_minus = copy(x_curr)
        x_minus[i] = x_minus[i] - Δ[i]

        x_plus = copy(x_curr)
        x_plus[i] = x_plus[i] + Δ[i]

        f_min = min(f(x_curr),min(f(x_minus),f(x_plus)))

        if f(x_minus) == f_min
            x_curr = x_minus
            success=true
        elseif f(x_plus) == f_min
            x_curr = x_plus
            success=true
        end
    end

    return x_curr, success
end

point0 = [-1.2,1]
Δ = [1.0,1.0]
τ = [0.000000001,0.000000001]
rosenbrock(x; a=1.0, b=100) = (a-x[1])^2 + b*(x[2]-x[1]^2)^2
result = direct_search_jeeve_hooke(rosenbrock,point0,Δ,τ)
println(result)
