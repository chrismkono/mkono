%   ---------------------------------------------------------------
%   Function Name:  expmodel

function expmodel(experimodel)

fid=fopen('exp_model.xls','w+');
for i=1:size(experimodel,1)
    fprintf(fid,'%d',i);
    for j=1:size(experimodel,2)
        fprintf(fid,'\t%17.10f',experimodel(i,j)); 
    end
    fprintf(fid,'\n');
end
fprintf(fid,'* * *\n');
fclose(fid);
fprintf(1,'\n The model results has been saved in ''exp_model.xls''.\n\n')

