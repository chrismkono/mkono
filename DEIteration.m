%   ---------------------------------------------------------------
%   Function Name:  DE Iteration

function DEPopulation = DEIteration(xin,pop,fit,seed,CR,F)
seed = ran2(-1e6*seed);
% DEPop - Replicate the population for operation in the function
DEPop = pop;

% ProbSize - gives the size of the problem
% ProbDmension - gives the dimension of the problem i.e. 3
[ProbSize,ProbDimension] = size(xin);

% PopSize - the number of individuals in the population
% PopDimension - the number of elements in the individual 
[PopSize,PopDimension] = size(DEPop);

% BestValue - the best value in the population
% BestPosition - the position of the best value in the population
[BestValue,BestPosition] = max(fit);

%Preallocate TrialPopulation
%TrialPopulation = zeros(1,PopSize);

% Forward Transformation - Convert the population to real domain
DERealPopulation = DEForwardTransformation(DEPop,PopSize,PopDimension);
DEPop = DERealPopulation;

for Iter = 1:PopSize
    if (Iter == BestPosition)
        TrialPopulation(Iter,:) = DEPop(Iter,:);
    else
        DEIndividual = DERoutine(DEPop,PopSize,PopDimension,BestPosition,Iter,seed,CR,F);
        TrialPopulation(Iter,:) = DEIndividual;
    end
end

% Backward Transformation - Convert the population back to real domain
DEDiscretePopulation = DEBackwardTransformation(TrialPopulation,PopSize,PopDimension);
DEPopulation = DEDiscretePopulation;

% Check for feasibility. If the bounds are violated, then move the
% offending elements to the boundary values.
for PopIter = 1:PopSize
    for IndiIter = 1:PopDimension
        if (DEPopulation(PopIter,IndiIter) < 1)
            seed = ran2(-1e6*seed);
            DEPopulation(PopIter,IndiIter) = ceil(ProbDimension*ran2(-1e6*seed));
        elseif (DEPopulation(PopIter,IndiIter) > ProbDimension)
            seed = ran2(-1e6*seed);
            DEPopulation(PopIter,IndiIter) = ceil(ProbDimension*ran2(-1e6*seed));
        end
    end
end