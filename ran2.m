%   ---------------------------------------------------------------
%   Function Name:  ran 2 - Random Number Generator


function ran2=ran2(idum)
IM1=2147483563;
IM2=2147483399;
AM=1./IM1;
IMM1=IM1-1; 

IA1=40014;
IA2=40692;
IQ1=53668;
IQ2=52774;
IR1=12211;
IR2=3791; 
NTAB=32;
NDIV=1+fix(IMM1/NTAB);
EPS=1.2e-7;
RNMX=1.-EPS;
idum2=123456789; 
iv=NTAB*0; 
iy=0;

if idum<=0
    idum=max(-idum,1);
    idum2=idum;
    for j=NTAB+8:-1:1
        k=fix(idum/IQ1);
        idum=IA1*(idum-k*IQ1)-k*IR1;
        if idum<0;idum=idum+IM1;
        end
        if j<=NTAB;iv(j)=idum;
        end
    end
    iy=iv(1);
end

k=fix(idum/IQ1);
idum=IA1*(idum-k*IQ1)-k*IR1;

if idum<0;idum=idum+IM1;end
k=fix(idum2/IQ2);
idum2=IA2*(idum2-k*IQ2)-k*IR2;

if idum2<0;idum2=idum2+IM2;end
j=1+fix(iy/NDIV);
iy=iv(j)-idum2;
iv(j)=idum;

if iy<1;iy=iy+IMM1;end

ran2=min(AM*iy,RNMX);
