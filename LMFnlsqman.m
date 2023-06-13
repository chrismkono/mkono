%   ---------------------------------------------------------------
%   Function Name:  LMFnlsqmanman

function [xf, SS, cnt, res, XY] = LMFnlsqman(varargin) 
echo off
if nargin==0 && nargout==0, help LMFnlsqman, return, end

if nargin==0 || (nargin==1 && strcmpi('default',varargin(1)))
    xf.Display  = 0;
    xf.Jacobian = 'finjac';
    xf.MaxIter  = 100;
    xf.ScaleD   = [];
    xf.FunTol   = 1e-7;
    xf.XTol     = 1e-7;
    xf.Printf   = 'printit';
    xf.Trace    = 0;
    xf.Lambda   = 0;
    return

elseif isstruct(varargin{1}) 
   if ~isfield(varargin{1},'Jacobian')
        error('Options Structure not Correct for LMFnlsqman.')
    end
    xf=varargin{1};         
    for i=2:2:nargin-1
        name=varargin{i}; 
        if ~ischar(name)
            error('Parameter Names Must be Strings.') 
        end
        name=lower(name(isletter(name))); 
        value=varargin{i+1}; 
        if strncmp(name,'d',1), xf.Display = value;
        elseif strncmp(name,'f',1), xf.FunTol = value(1);
        elseif strncmp(name,'x',1), xf.XTol = value(1);
        elseif strncmp(name,'j',1), xf.Jacobian = value;
        elseif strncmp(name,'m',1), xf.MaxIter = value(1);
        elseif strncmp(name,'s',1), xf.ScaleD = value;
        elseif strncmp(name,'p',1), xf.Printf = value;
        elseif strncmp(name,'t',1), xf.Trace = value;
        elseif strncmp(name,'l',1), xf.Lambda = value;
        else disp(['Unknown Parameter Name --> ' name]) 
        end
    end
    return

elseif ischar(varargin{1}) 
    Pnames=char('display','funtol','xtol','jacobian','maxiter','scaled',... 
        'printf','trace','lambda');
    if strncmpi(varargin{1},Pnames,length(varargin{1})) 
        xf=LMFnlsqman('default'); 
        xf=LMFnlsqman(xf,varargin{:});
        return
    end
end

FUN=varargin{1}; 
if ~(isvarname(FUN) || isa(FUN,'function_handle'))
   error('FUN Must be a Function Handle, or M-file Name.')
end
xc=varargin{2};             
if ~exist('options','var')
    options = LMFnlsqman('default');
end
if nargin>2                 
    if isstruct(varargin{3})
        options=varargin{3};
    else
        for i=3:2:size(varargin,2)-1
            options=LMFnlsqman(options,varargin{i},varargin{i+1});
        end
    end
else
    if ~exist('options','var')
        options = LMFnlsqman('default');
    end
end

x = xc(:);
n = length(x);
epsx = options.XTol(:);
le = length(epsx);
if le==1
    epsx=epsx*ones(n,1);
else
    error(['Dimensions of vector epsx ',num2str(le),'~=',num2str(lx)]); 
end
epsf  = options.FunTol(:);
ipr   = options.Display;
JAC   = options.Jacobian;
maxit = options.MaxIter;   
printf= options.Printf;

r = feval(FUN,x);
[A,v] = getAv(FUN,JAC,x,r,epsx);

SS = r'*r;
res= 1;
cnt=0;
trcXY = options.Trace;
if trcXY
    XY = zeros(n,maxit);
    XY(:,1) = x;
else
    XY = []; 
end

D = options.ScaleD(:); 
if isempty(D)   
    D=diag(A); 
else
    ld=length(D);
    if ld==1
        D=abs(D)*ones(n,1); 
    elseif ld~=n
        error(['wrong number of scales D, lD = ',num2str(ld)])
    end
end
D(D<=0)=1;
T = sqrt(D);
Rlo=0.25;           Rhi=0.75;
l=options.Lambda;   lc=1;
dx = zeros(n,1);
cnt = 0;

while 1 
    feval(printf,ipr,cnt,res,SS,x,dx,l,lc)
    cnt = cnt+1;
    if trcXY, XY(:,cnt+1)=x; end
 
    d = diag(A);
    s = zeros(n,1);

    while 1 
        while 1
            UA = triu(A,1);
            A = UA'+UA+diag(d+l*D);
            [U,p] = chol(A);        %
            %~~~~~~~~~~~~~~~
            if p==0, break, end
            l = 2*l;
            if l==0, l=1; end
        end
        dx = U\(U'\v);
        vw = dx'*v;
        fin = -1;
        if vw<=0, break, end       
        
        for i=1:n
            z = d(i)*dx(i);
            if i>1, z=A(i,1:i-1)*dx(1:i-1)+z; end 
            if i<n, z=A(i+1:n,i)'*dx(i+1:n)+z; end 
            s(i) = 2*v(i)-z;
        end
        dq = s'*dx;
        s  = x-dx;
        rd = feval(FUN,s);
             
        res = res+1;
        SSP = rd'*rd;
        dS = SS-SSP;
        fin = 1;
        if all((abs(dx)-epsx)<=0) || res>=maxit || abs(dS)<=epsf
            break                  
        end
        fin=0;
        if dS>=Rlo*dq, break, end
        A = U;
        y = .5;
        z = 2*vw-dS;
        if z>0, y=vw/z; end
        if y>.5, y=.5; end
        if y<.1, y=.1; end
        if l==0
            y = 2*y;
            for i = 1:n
                A(i,i) = 1/A(i,i);
            end
            for i = 2:n
                ii = i-1;
                for j= 1:ii
                    A(j,i) = -A(j,j:ii)*A(j:ii,i).*A(i,i);
                end
            end
            for i = 1:n
                for j= i:n
                    A(i,j) = abs(A(i,j:n)*A(j,j:n)');
                end
            end
            l  = 0;
            tr = diag(A)'*D;
            for i = 1:n
                z = A(1:i,i)'*T(1:i)+z;
                if i<n
                    ii = i+1;
                    z  = A(i,ii:n)*T(ii:n)+z;
                end
                z = z*T(i);
                if z>l, l=z; end
            end
            if tr<l, l=tr; end
            l  = 1/l;
            lc = l;
        end
        l = l/y;
        if dS>0, dS=-1e300; break, end
    end 
    
    if fin, break, end
    if dS>Rhi*dq
        l=l/2;
        if l<lc, l=0; end
    end
    SS=SSP;  x=s;  r=rd;
    [A,v] = getAv(FUN,JAC,x,r,epsx);

end 

if fin>0
    if dS>0
        SS = SSP;
        x = s; 
    end
end
if ipr~=0
    disp(' ');
    feval(printf,sign(ipr),cnt,res,SS,x,dx,l,lc) 
end
xf = x;
if trcXY, XY(:,cnt+2)=x; end
XY(:,cnt+3:end) = [];
if res>=maxit, cnt=-maxit; end
return 
