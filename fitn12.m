%   ---------------------------------------------------------------
%   Function Name:  fitn12

function para=fitn12(xin,len,pop,qwe,y,qwe1,metofcal,ninp); 

global IG
des_var=pop;
len_fun=8;

[newsizex,newsizey]=size(pop);
npop=newsizex;
xreg=xin;
[cx,xc]=size(xin);
for iis=1:npop
    for ksk=1:2:len(iis)-1
        if pop(iis,ksk)>pop(iis,ksk+1)
            sis=pop(iis,ksk+1);
            pop(iis,ksk+1)=pop(iis,ksk);
            pop(iis,ksk)=sis;
        end
    end
end
xpop=pop;
for i=1:npop
    xin=xreg;
    n=len(i);
    obj_Layer(i)=log(n)/log(2);
    if n==2
        neron(1,1:2)=pop(i,1:2);
    else
        layer=log(n)/log(2);
        lay=1;
        neron(lay,:)=pop(i,:);
        while lay<=(layer-1)
            qw=1;
            for oiu=1:2:len(i)-1
                neron(lay+1,qw)=neron(lay,oiu)*((ninp+1)^(2^(lay- 1)))+neron(lay,oiu+1);
                qw=qw+1; 
            end
            for ksk=1:2:len(i)-1
                if neron(lay+1,ksk)>neron(lay+1,ksk+1)
                sist=neron(lay+1,ksk+1); 
                neron(lay+1,ksk+1)=neron(lay+1,ksk); 
                neron(lay+1,ksk)=sist;
                end
            end
            lay=lay+1;
        end
    end
    [azz,bzz]=size(neron);
    neuron_countor=0;
    counting=0;
    for j=1:azz
        jj=1;
        for uu=1:2:len(i)
            if j==1
                u=neron(j,uu);
                v=neron(j,uu+1);
            else
                u=uu;
                v=uu+1; 
            end
            if neron(j,uu)==neron(j,uu+1)
                xnew(:,jj)=xin(:,u);
                if j~=1
                    counting=counting+1;
                end
            else
                neuron_countor=neuron_countor+1; 
                for w=1:qwe
                    A(w,:)=[1,xin(w,u),xin(w,v),xin(w,u)^2,xin(w,v)^2,xin(w,u)*xin(w,v)]; 
                end
                if metofcal==1
                    a=pinv(A'*A)*A'*(y(1:qwe))';
                    a=a';
                elseif metofcal==2
                    [AAa,BBb,CCc]=svds(A);
                    [hH,gG]=size(BBb);
                    for iII=1:gG
                        if BBb(iII,iII)>1e-3
                            BBb(iII,iII)=1/BBb(iII,iII);
                        else
                            BBb(iII,iII)=0;
                        end
                    end
                    a=CCc*BBb*AAa'*(y(1:qwe))';
                    a=a';
                elseif metofcal==3
                    a=pinv(A'*A)*A'*(y(1:qwe))';
                    a=a';
                    [x]=run_lm(xin,qwe,A,u,v,a);
                    a=x;
                elseif metofcal==4
                    a=[IG]';
                    [x]=run_lm(xin,qwe,A,u,v,a);
                    a=x; 
                end
                clear CCc BBb AAa
                for w=1:qwe1
                    xnew(w,jj)=(a(1,1)+a(1,2)*xin(w,u)+a(1,3)*xin(w,v)+...
                        +a(1,4)*xin(w,u)^2+a(1,5)*xin(w,v)^2+a(1,6)*xin(w,u)*xin(w,v)); 
                end
            end
            jj=jj+1; 
        end
        [msize,nsize]=size(neron);
        countour=0;
        nsize=nsize/2;
        pop_man=pop(i,:);
        for m_ner=2:msize
            nfind=0;
            for n_ner=1:nsize
                if (m_ner==1)&(pop_man(m_ner-1,2*n_ner-1)~=pop_man(m_ner-1,2*n_ner))
                    if (n_ner~=1) 
                        [mfind,nfind]=size(find(neron(m_ner,1:n_ner-1)==neron(m_ner,n_ner)));
                    end
                    if nfind==0
                        [m_count,n_count]=size(find(neron(m_ner,:)==neron(m_ner,n_ner))); 
                        if n_count>1
                            countour=countour+n_count-1;
                        end
                    end
                end
            end
            nsize=nsize/2;
        end
        neuron_fun(i)=neuron_countor-countour; 
        len(i)=len(i)/2;
        clear xin
        xin=xnew;
        clear xnew
    end
    clear neron a
    err(i)=0;
    for we=1:qwe
        err(i)=err(i)+(xin(we)-y(we))^2; 
    end

    err_predict(i)=0;
    for wep=qwe+1:qwe1
        err_predict(i)=err_predict(i)+(xin(wep)-y(wep))^2; 
    end

    yout(:,i)=xin;
end

obj_fun(1,:)=err./qwe; 
obj_fun(2,:)=err_predict./(qwe1-qwe);

para=[des_var obj_fun'];