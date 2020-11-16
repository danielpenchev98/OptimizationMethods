using DataStructures

struct Point2D
    x::Float64
    y::Float64
end

function checkIfCounterClockwiseTurn2D(a,b,c)
    return a.x*b.y + b.x*c.y + a.y*c.x - c.x*b.y - c.y*a.x - b.x*a.y > 0
end

function calculateAngleInPolarCoord(a,b)
    Δx = b.x - a.x
    Δy = b.y - a.y
    return atan(Δy,Δx) * 180.0 / π
end

#X
function convexHull(points)

    start = points[1]
    for i in 2:size(points,1)
        if start.y > points[i].y
            start = points[i]
        end
    end

    sort!(points, by = p -> calculateAngleInPolarCoord(start,p))

    path = Stack{Point2D}()
    push!(path,start)
    push!(path,points[2])

    for i in 3:size(points,1)
        curr = pop!(path)

        if !checkIfCounterClockwiseTurn2D(first(path),curr,points[i])
            push!(path,points[i])
            continue
        end

        push!(path,curr)
        push!(path,points[i])

    end

    final = pop!(path)
    if checkIfCounterClockwiseTurn2D(first(path),final,start)
        push!(path,final)
    end
    while !isempty(path)
        println(pop!(path))
    end
end

points = Array{Point2D}(undef, 6)
points[1] = Point2D(0.0,0.0)
points[2] = Point2D(5.0,1.0)
points[3] = Point2D(0.0,5.0)
points[4] = Point2D(-5.0,3.0)
points[5] = Point2D(-1.0,1.0)
points[6] = Point2D(1.0,0.5)

convexHull(points)
