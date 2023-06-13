%   ---------------------------------------------------------------
%   Function Name:  muta


function mutpop=muta(pop,len,npop,ninp,pmm,maxlen,sortpop,prob,seed) 
seed=ran2(-1e6*seed);
mutpop=sortpop;
nmut=ceil(pmm*npop);
for j=1:nmut
    seed=ran2(-1e6*seed);
    mupop(j)=ceil(npop*ran2(-1e6*seed)); 
end
for zx=1:nmut
    for ia=1:len(zx)
        seed=ran2(-1e6*seed);
        poi=ceil(ninp*ran2(-1e6*seed));
        mutpop(zx,ia)=poi;
    end
end