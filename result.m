%   ---------------------------------------------------------------
%   Function Name:  results

function result(bestfit)

fid=fopen('bestfitness.xls','w+');
for i=1:size(bestfit,2)
    fprintf(fid,'%d',i);
        fprintf(fid,'\t%17.10f',bestfit(i));
    fprintf(fid,'\n');
end
fprintf(fid,'* * *\n');
fclose(fid);
fprintf(1,'\n The bestfitness has been saved in ''bestfitness.xls''.\n\n')