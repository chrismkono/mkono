%   ---------------------------------------------------------------
%   Function Name:  propop

function [sortpop,prob]=propop(npop,fit,sortfit,pop,maxlen) 
i=1;
while i<=npop
    [oo,pp]=find(fit==sortfit(i));
    [ii,jj]=size(pp);
    for b=1:jj
        sortpop(i,1:maxlen)=pop(pp(b),1:maxlen);
        i=i+1; 
    end
end
tofit=sum(fit);
prob=sortfit/tofit;

prob=cumsum(prob);