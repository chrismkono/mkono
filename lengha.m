%   ---------------------------------------------------------------
%   Function Name:  lengha

function len=lengha(pop,npop)
for oi=1:npop
    [om,on]=size(find(pop(oi,:)));
    len(oi)=on;
end