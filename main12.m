%   ---------------------------------------------------------------
%   Function Name:  main12 - Main Function

function [y,bestfit,final]=main12(data,N_pop,N_HL,pmm,N_Generation,N_Prediction,metofcal,first_seed,CR,F)

seed=ran2(first_seed); 
[Rinputs,Cinputs]=size(data); 

for i=1:Cinputs-2
    xin(:,i)=data(:,i+1);
end
y=data(:,Cinputs);
y=y';
[qwe1,mnb]=size(xin);

ninp=mnb;

npop=N_pop;
maxlen=2^(N_HL+1);
niter=N_Generation;
pre=N_Prediction;
qwe=qwe1-pre;
nlay=log(maxlen)/log(2);

seed=ran2(-1e6*first_seed);

for i=1:npop
    seed=ran2(-1e6*seed);
    len(i)=maxlen*ran2(-1e6*seed);
end

k=maxlen;

while k>1
    b=find(len<k&len>k/2);
    len(b)=k;
    k=k/2;
    clear b
end

bb=find(len<=2);
len(bb)=2;
seed=ran2(-1e6*seed);

for i=1:npop
    seed=ran2(-1e6*seed);
    for j=1:len(i)
        seed=ran2(-1e6*seed);
        s(j)=ceil(ninp*ran2(-1e6*seed)); 
    end
    
    bB=find(s==s(1));
    [m,n]=size(s);
    
    [mm,nn]=size(bB);
    if n==nn
        if s(1)==ninp
            s(1)=1;
        else
            s(1)=s(1)+1;
        end
    end
    pop(i,1:n)=s;
    for k=1:2:len(i)-1
        if pop(i,k)>pop(i,k+1)
            ss=pop(i,k+1);
            pop(i,k+1)=pop(i,k);
            pop(i,k)=ss;
        end
    end
end

inipop=pop;
xreg=xin;
seed=ran2(-1e6*seed);

for iter=1:niter    
    seed=ran2(-1e6*seed);
    xin=xreg; 

    [yout,xpop,fit,sortfit,averfit]=fitn2(xin,len,npop,pop,qwe,y,qwe1,metofcal,ninp);

    averfitt(iter)=averfit;
    pop=xpop;
    xin=xreg;
    [uu,tt]=max(fit);

    bestfit(iter)=max(fit);

    fina(iter,1:maxlen)=pop(tt,1:maxlen); 
    [sortpop,prob]=propop(npop,fit,sortfit,pop,maxlen); 
    seed=ran2(-1e6*seed);

    repop=repro(npop,sortpop,prob,seed);
    [sortpop,prob]=propop(npop,fit,sortfit,repop,maxlen); 

    xin=xreg;
    seed=ran2(-1e6*seed); 

    DEPopulation = DEIteration(xin,pop,fit,seed,CR,F);
    pop = DEPopulation;
    
    len=lengha(pop,npop);

    [yout,xpop,fit,sortfit,averfit]=fitn2(xin,len,npop,pop,qwe,y,qwe1,metofcal,ninp);

    [sortpop,prob]=propop(npop,fit,sortfit,pop,maxlen);

    mutpop=muta(pop,len,npop,ninp,pmm,maxlen,sortpop,prob,seed);  

    mutpop(npop,1:maxlen)=fina(iter,1:maxlen); 
    pop=mutpop;
    len=lengha(pop,npop);
end

xin=xreg; 
[yout,xpop,fit,sortfit,averfit]=fitn2(xin,len,npop,pop,qwe,y,qwe1,metofcal,ninp);
[uu,tt]=max(fit);
averfitt(iter+1)=averfit;
ynout(:,1)=yout(:,tt);

pop=xpop;
bestfit(niter+1)=max(fit); 
finalpop(niter+1,1:maxlen)=pop(tt,1:maxlen); 
xin=xreg;

clc
fprintf('\n\n Best Population for this Experimental Results is:\n '); 
finalpop(niter+1,:)
final=finalpop(niter+1,:);
bestpop(finalpop(niter+1,:));
fprintf('\n\n Best Fitness equal to :');
bestfit(niter+1)
hor=1:niter+1;
result(bestfit)

experimodel(:,1)=y'; 
experimodel(:,2)=ynout; 
expmodel(experimodel) 
lil=min(y):(max(y)-min(y))/100:max(y); 
[poi,iop]=size(lil); 
prt=qwe*ones(1,iop);

fprintf('\n\n Neural Network Multiplyer for Final population:\n '); 
outpop=finalpop(niter+1,:); 

format long

global data
global N_Prediction

[Rinputs,Cinputs]=size(data);
for i=1:Cinputs-2
    xin(:,i)=data(:,i+1);
end
y=data(:,Cinputs);
y=y';
[qwe1,mnb]=size(xin);
Xin_col=mnb;
N_trn=qwe1-N_Prediction;
Best_pop=[
    final
    ]
yout=multip(xin,y,Xin_col,Best_pop,N_trn,metofcal);
yout=yout'
Re=(abs(y-yout))./y 

MRe_Train=sum(Re(1,1:10))/10
MRe_Predict=sum(Re(1,11:15))/5 

para=fitn12(xin,8,Best_pop,N_trn,y,qwe1,1,mnb); 
plot(1:qwe1,y,'-o',1:qwe1,yout','-*')

toc
