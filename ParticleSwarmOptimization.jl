using Distributions

function particleSwarmOptimization()
    #initialize the population
end

function generatePopulation(populationSize, dimensions,dimensionCostraints)
    population = [zeros(dimensions) for _ in 1:populationSize]
    for i in 1:populationSize
        for j in 1:dimension
            population[i][j] = rand(Uniform(dimensionConstraints[j][1]:dimensionConstraints[j][2]))
        end
    end

    return generatePopulation
end


 dimensionConstraints = [ [-5.0,5.0], [-5.0,-5.0] ]

 population = generatePopulation(10,size(dimensionConstraints)[1],dimensionConstraints)

 for i in 1:10
     println(population[i])
 end
