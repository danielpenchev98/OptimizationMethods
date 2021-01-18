using Random
using Distributions
using StatsBase
using Plots
function geneticAlgorithm(generationSize,itemValues,itemWeights,maxWeight,threshold;generationLimit=60)
    chromosomeSize = size(itemValues)[1]
    currGeneration = [zeros(chromosomeSize) for _ in 1:generationSize]
    fitness = zeros(generationSize)

    maxFitness = 0
    solution = zeros(chromosomeSize)

    for i in 1:generationSize
        for j in 1:chromosomeSize
            currGeneration[i][j]=rand(0:1)
        end
        fitness[i] = getFitness(currGeneration[i],itemValues,itemWeights,maxWeight)
        if maxFitness < fitness[i]
            maxFitness = fitness[i]
            solution = currGeneration[i]
        end
    end

    recordMaxFitnees = zeros(generationLimit)
    recordMaxFitnees[1]=maxFitness


    generationNum = 1
    while maxFitness < threshold && generationNum <= generationLimit
        maxFitness=0
        currGeneration = createNewGeneration(currGeneration,fitness,chromosomeSize)
        for i in 1:generationSize
            fitness[i] = getFitness(currGeneration[i],itemValues,itemWeights,maxWeight)
            if maxFitness < fitness[i]
                maxFitness = fitness[i]
                solution = currGeneration[i]
            end
        end
        recordMaxFitnees[generationNum]=maxFitness
        generationNum+=1
    end
    println(solution)

    plot(1:size(recordMaxFitnees)[1],recordMaxFitnees)
end

function getFitness(chromosome,itemValues,itemWeights,maxWeight)
    totalSum, totalWeight = 0, 0
    for i in 1:size(chromosome)[1]
        if chromosome[i] == 1.0
            totalSum+=itemValues[i]
            totalWeight+=itemWeights[i]
        end
    end

    if(totalWeight>maxWeight)
        totalSum=0
    end
    return totalSum
end

function createNewGeneration(currGeneration,fitness,sizeOfChromosome; numCompetitors=30)
    println(size(currGeneration)[1])
    sizeOfGeneration = size(currGeneration)[1]

    winners = [zeros(sizeOfChromosome) for _ in 1:sizeOfGeneration]

    #schedule tournaments
    tournaments = [sample(1:sizeOfGeneration, numCompetitors ,replace = false) for _ in 1:sizeOfGeneration]

    #play tournaments
    for i in 1:sizeOfGeneration
        winnerIndex = getWinnerFromTournament(tournaments[i],fitness)
        winners[i] = currGeneration[winnerIndex]
    end

    #breed new generation
    newGeneration = [zeros(sizeOfChromosome) for _ in 1:sizeOfGeneration]

    i = 1
    while i < sizeOfGeneration
        newGeneration[i], newGeneration[i+1] = getChildren(winners[i],winners[i+1])
        i+=2
    end

    #if population size is odd number, than the last winner will mutate and be added to next generation
    if i==sizeOfGeneration
        mutate(winners[i])
        newGeneration[i] = winners[i]
    end

    return newGeneration
end

function getWinnerFromTournament(choosen,fitness)
    best = choosen[1]
    for i in choosen
        if fitness[best] < fitness[i]
            best = i
        end
    end
    return best
end

function getChildren(firstParent,secondParent)
    firstChild, secondChild = crossover(firstParent,secondParent);
    mutate(firstChild)
    mutate(secondChild)
    return firstChild, secondChild
end

function crossover(firstParent,secondParent)
    len = size(firstParent)[1]
    firstChild = zeros(len)
    secondChild = zeros(len)
    for i in 1:len#in Int(len/2):len
        if i<=floor(Int,len/2)
            firstChild[i],secondChild[i] = firstParent[i], secondParent[i]
        else
            firstChild[i],secondChild[i] = secondParent[i], firstParent[i]
        end
    end
    return firstChild, secondChild
end

function mutate(chromosome; mutationProb=0.05)
    distribution = Binomial(1,mutationProb)

    for i in 1:size(chromosome)[1]
        shouldMutate = rand(distribution,1)[1]
        if shouldMutate == 1
            chromosome[i] = (chromosome[i]+1)%2
        end
    end
end


#mutate(ones(100),0.05)

#geneticAlgorithm(3,10)
#f1 = [0.0,0.0,1.0,0.0]
#f2 = [1.0,0.0,0.0,1.0]

#println(crossover(f1,f2))

#itemValues = [12.0,0.2,0.4,0.5]
#itemWeights = [2.3,1.0,0.3,0.1]
#maxWeight = 3.0

#items = [1.0,1.0,0.0,1.0]
#println(getFitness(items,itemValues,itemWeights,maxWeight))


#fitness = [12.0, 1.5, 134.6 , 1.7]
#fitness = [-1000.7,1.5,12.0,134.6]
#println(getWinnerFromTournament([1,2,3,4],fitness))

#firstParent = [1.0,0.0,1.0,1.0,0.0,1.0]
#secondParent = [0.0,1.0,0.0,0.0,1.0,0.0]
#println(getChildren(firstParent,secondParent))

#f1 = [1.0,0.0,1.0,1.0,0.0,0.0]
#f2 = [0.0,1.0,0.0,1.0,1.0,1.0]
#f3 = [1.0,0.0,0.0,0.0,1.0,0.0]
#f4 = [0.0,1.0,1.0,0.0,0.0,0.0]
#f5 = [0.0,1.0,0.0,0.0,0.0,1.0]

#itemValues = [5.0,2.0,6.0,1.0,2.0,3.0]
#itemWeights = [4.0,1.0,7.0,3.0,2.0,5.0]

chromosomes = [zeros(40) for _ in 1:100]
itemValues = rand(0:40,40)
itemWeights = rand(0:40,40)

println(itemValues)

maxWeight = 500

geneticAlgorithm(100,itemValues,itemWeights,maxWeight,800)
