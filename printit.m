%   ---------------------------------------------------------------
%   Function Name:  printit

function printit(ipr,cnt,varargin)

echo off
if ipr~=0 && rem(cnt,ipr)==0
    lv  = length(varargin);     
    if lv>4 && ipr<0
        lv=lv-2;                
    end
    n = min(75-26*(ipr<0),13*lv-3);
    hlin = @(n) fprintf(['\n',repmat('*',1,n),'\n']);
    if cnt==0                   
        fulh ={' itr',' nfJ',' sum(r^2)', ' x',...
            ' dx',' lambda',' lambda_c'...
            }; 
        hlin(n);
        fprintf('%s',fulh{1:lv+1});        
        hlin(n);
    end
    lx  = min(4,lv);                
    xdx = [varargin{3:lx}];
    var = [varargin{1:2},xdx(1,:)]; 
    if lv>4 && ipr>0
        var = [var varargin{5:6}];  
    end
    fprintf(['%4.0f %4.0f ' repmat('%12.4e ',1,lv-1),'\n'], cnt,var); 
    fprintf([blanks(23),repmat('%12.4e ',1,lx-2),'\n'], xdx(2:end,:)'); 
end