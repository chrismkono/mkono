%   ---------------------------------------------------------------
%   Function Name:  DE Routine

function DEIndividual = DERoutine(DEPop,PopSize,PopDimension,BestPosition,Iter,seed,CR,F)
seed = ran2(-1e6*seed);

% NumOfVectors  -  number of individuals selected for crossover
NumOfVectors = 2;

% VectorIndex - array to hold the number of individuals selected for
% crossover + best individual
VectorIndex = [0 0];

% Obtain two random index for the DE crossover operation pairwise different
% from each other, the best solution and the current index.
for Index = 1:NumOfVectors
    seed=ran2(-1e6*seed);
    RandomIndex = ceil(PopSize*ran2(-1e6*seed));
    if Index == 1
        while((RandomIndex == Iter) || (RandomIndex == BestPosition))
            seed=ran2(-1e6*seed);
            RandomIndex = ceil(PopSize*ran2(-1e6*seed));
        end
    elseif Index > 1
        while(RandomIndex == Iter || RandomIndex == BestPosition || RandomIndex == VectorIndex(1))
            seed=ran2(-1e6*seed);
            RandomIndex = ceil(PopSize*ran2(-1e6*seed));
        end
    end    
    VectorIndex(Index) = RandomIndex;
end

% Obtain a unique start position
seed = ran2(-1e6*seed);

% Set and staart the counter for iterating the individual. As the first
% value is always changed, Counter is set to 2.
Counter = 2;
StartPosition = ceil(PopDimension*ran2(-1e6*seed));

% Copy the index individual to temporary array for evalaution
DEIndividual = DEPop(Iter,:);

% Compute the first element in the new trial individual
DEIndividual(StartPosition) = DEPop(BestPosition,StartPosition) + (F * (DEPop(VectorIndex(1),StartPosition) - DEPop(VectorIndex(2),StartPosition)));

% Check is the StratPosition is not at the end of the Individual
if(mod(StartPosition,PopDimension) == 0)
    StartPosition = 1;
else
    StartPosition = StartPosition + 1;
end

seed = ran2(-1e6*seed);

% Iterate through the individual while the CR condition is specified.
while ((ran2(-1e6*seed) < CR) && (Counter <= PopDimension))
    DEIndividual(StartPosition) = DEPop(BestPosition,StartPosition)+ (F * (DEPop(VectorIndex(1),StartPosition) - DEPop(VectorIndex(2),StartPosition)));
    
    % Increment the start position
    if(mod(StartPosition,PopDimension) == 0)
        StartPosition = 1;
    else
        StartPosition = StartPosition + 1;
    end
    
    % Increment the counter variable
    Counter = Counter + 1;
end
