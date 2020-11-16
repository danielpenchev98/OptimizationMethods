using Distributions

function particleSwarmOptimization(f,dimensions,dimensionsConstraints,populationSize; startWeight=0.9,minWeight=0.4,c1=2.0,c2=2.0,maxIter=100, maxVelocity=20.0)
    #initialize the population
    particlePositions = generatePopulation(populationSize,size(dimensionsConstraints)[1],dimensionsConstraints)

    #initalize the personal best
    personalBestPositions = deepcopy(particlePositions)

    #initialize global best
    globalBestPosition = particlePositions[1]
    globalBestPosFitness = f(globalBestPosition)
    for i in 2:populationSize
        currValue = f(particlePositions[i])
        if globalBestPosFitness > currValue
            globalBestPosFitness = currValue
            globalBestPosition = particlePositions[i]
        end
    end

    #initialize velocities
    velocities = initVelocities(startWeight,particlePositions,dimensions,personalBestPositions,globalBestPosition,c1,c2)

    currWeight = startWeight
    weightUpdateVal = (startWeight - minWeight)/maxIter

    currIter = 1
    while currIter <= maxIter
        for i in 1:populationSize
            #update personal best and global best
            currFitness = f(particlePositions[i])
            if f(personalBestPositions[i]) > currFitness
                personalBestPositions[i] = particlePositions[i]
            end

            if f(globalBestPosition) > currFitness
                globalBestPosition = particlePositions[i]
            end
        end

        #update weight
        currWeight = currWeight - weightUpdateVal

        #update position of particle
        for i in 1:populationSize
            velocities[i] = getNextVelocity(currWeight,velocities[i],particlePositions[i],personalBestPositions[i],globalBestPosition,c1,c2)
            particlePositions[i] = getNextPosition(particlePositions[i] ,velocities[i])
        end

        currIter = currIter+1
    end

    return globalBestPosition
end

function initVelocities(weight,particlePositions,dimensions,personalBestPos,globalBestPos,c1,c2)
    velocities = [zeros(dimensions) for _ in 1:size(particlePositions)[1]]
    for i in 1:size(particlePositions)[1]
        velocities[i]=getNextVelocity(weight,velocities[i],particlePositions[i],particlePositions[i],globalBestPos,c1,c2)
    end
    return velocities
end

function getNextVelocity(weight,prevVelocity,position,personalBest,globalBest,c1,c2)
    r1 = rand(Uniform(0.0,0.1),1)[1]
    r2 = rand(Uniform(0.0,0.1),1)[1]
    inertia = weight*prevVelocity
    cognitiveComponent = c1 * r1 * (personalBest - position)
    socialComponent = c2 * r2 * (globalBest - position)
    return inertia +  cognitiveComponent + socialComponent
end

function getNextPosition(currPosition,newVelocity)
    return currPosition + newVelocity
end

function generatePopulation(populationSize, dimensions,dimensionCostraints)
    population = [zeros(dimensions) for _ in 1:populationSize]
    for i in 1:populationSize
        for j in 1:dimensions
            population[i][j] = rand(Uniform(dimensionConstraints[j][1],dimensionConstraints[j][2]),1)[1]
        end
    end

    return population
end


branin(x; a=1, b=5.1/(4π^2), c=5/π, r=6, s=10, t=1/(8π)) = a*(x[2]-b*x[1]^2+c*x[1]-r)^2 + s*(1-t)*cos(x[1]) + s
dimensionConstraints = [ [-10.0,20.0], [-10.0,20.0] ]

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

bestSolution = particleSwarmOptimization(branin,size(dimensionConstraints)[1],dimensionConstraints,40)

println("Best point so far :")
println(bestSolution)
scatter!([bestSolution[1]],[bestSolution[2]],label="Best solution")
