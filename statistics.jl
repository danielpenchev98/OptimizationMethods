using Distributions
using PyPlot
using CSV
using DataFrames
using GLM
using HypothesisTests
using StatsBase
using StatPlots
pyplot()

age = rand(18:80,100);
wcc = round.(rand(Distributions.Normal(12,2),100),digits=1);
crp = round.(Int,rand(Distributions.Chisq(4),100)) .* 10;
treatment = rand(["A","B"],100);
result = rand(["Improved","Static","Worse"],100);

#StatsBase.describe(age)
StatsBase.summarystats(age)

data = DataFrame(Age = age, WCC = wcc, CRP = crp, Treatment = treatment,Result = result);
print(size(data))
print(head(data))

dataA = data[data[:Treatment].== "A", :]
dataB = data[data[:Treatment] .== "B", :]

by(data, :Treatment,df ->DataFrame(N=size(df,1)))

by(data, :Treatment,size)

by(data, :Treatment, cf->mean(cf.Age))

#group by Treatment
@df data density( :Age, group = (:Treatment,:Result), title = "Distribution of ages by group type",xlab="Age",ylab="Distribution", legend= :topright)

@df data StatPlots.boxplot(:Treatment,:WCC,lab="WCC", title="White cells by treatment group", xlab="Groups",ylab="WCC")

@df data corrplot([:Age :WCC :CRP], grid=false)

@df data cornerplot([:Age :WCC :CRP], grid=false,compact=true)


## Inferential statistics - tests

HypothesisTests.EqualVarianceTTest(dataA[:Age],dataB[:Age])

pvalue(HypothesisTests.EqualVarianceTTest(dataA[:WCC],dataB[:WCC]))

HypothesisTests.UnequalVarianceTTest(dataA[:CRP],dataB[:CRP])

#Linear Regression

GLM.fit(LinearModel,@formula(CRP~1),data)


GLM.fit(LinearModel,@formula(CRP~Age),data)


GLM.fit(LinearModel,@formula(CRP~Age+WCC),data)

#Independance of categorial varialbe test - chisq

by(dataA,:Result,df->DataFrame(N=size(df,1)))
by(dataB,:Result,df->DataFrame(N=size(df,1)))

observed = reshape([22,7,19,19,15,18],(2,3))
ChisqTest(observed)
