using Convex
function strong_backtracking_line_search(f, ∇f, x0, d; α=1, β=exp(-4),σ = 0.1)
    g0 = dot(∇f(x0),d)
    y_prev, αlow, αhigh, αprev  =  NaN, NaN, NaN, 0
    ∇f0 = ∇f(x0)
    while true
        y = f(x0+α*d)
        if y > f(x0) + β*α*g0 || (!isnan(y_prev) && y >= y_prev)
            αlow, αhigh = αprev, α
            break
        end
        g = dot(∇f(x0 + α*d),d)
        if abs(g) <= -σ*g0
            return α
        elseif g >= 0
            αhigh, αlow = αprev, α
            break
        end
        αprev, α, y_prev = α, 2*α, y
    end

    #zoom phase
    while true
        y_low = f(x0 + αlow*d)
        α = (αlow + αhigh)/2
        y = f(x0 + α*d)
        if y > f(x0) + β*α*g0 || y >= y_low
            αhigh = α
        else
            g = dot(∇f(x0 + α*d),d)
            if abs(g) <= -σ*g0
                return α
            elseif g*(αhigh - αlow) >= 0
                αhigh = αlow
            end
            αlow = α
        end
    end
end

#using ReverseDiff

#x0 = [-1,-1]
#direction = [1/sqrt(2),1/sqrt(2)]
#rosenbrock(x) = (1-x[1])^2 + 105*(x[2]-x[1]^2)^2
#∇f = x -> ReverseDiff.gradient(rosenbrock,x)
#println(strong_backtracking_line_search(rosenbrock,∇f,x0,direction,α=1,σ=0.1))
