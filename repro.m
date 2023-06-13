%   ---------------------------------------------------------------
%   Function Name:  repro

function repop=repro(npop,sortpop,prob,seed)
seed=ran2(-1e6*seed);
% randomly generate a sequence of numbers (size of population)
for i=1:npop
    seed=ran2(-1e6*seed);
    rNums(i)=seed;
end
% sort those numbers in ascending order.
rNums=sort(rNums);

fitIn=1;newIn=1;
while newIn<=npop
    if(rNums(newIn)<prob(fitIn)) 
        repop(newIn,:) = sortpop(fitIn,:); 
        newIn = newIn+1;
    else
        fitIn = fitIn + 1;
    end
end