%   ---------------------------------------------------------------
%   Function Name:  bestpop

function bestpop(finalpop)

fid=fopen('bestpop.txt','w+');
fprintf(fid,'\t%4.0f',finalpop);
fprintf(fid,'\n * * *\n');
fclose(fid);
fprintf(1,'\n The best population has been saved in ''bestpop.txt''.\n\n')
