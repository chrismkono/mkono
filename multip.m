%   ---------------------------------------------------------------
%   Function Name:  multip

function yout=multip(xin,y,Xin_col,Best_pop,N_trn,metofcal) 
global IG
n_neroun=1;
[Sr,Sc]=size(Best_pop);
len=Sc;
xreg=xin;
[cx,xc]=size(xin);
qwe1=cx;
qwe=N_trn;
npop=Sr;
for iis=1:npop
    for ksk=1:2:len(iis)-1
        if Best_pop(iis,ksk)>Best_pop(iis,ksk+1)
        sis=Best_pop(iis,ksk+1); 
        Best_pop(iis,ksk+1)=Best_pop(iis,ksk); 
        Best_pop(iis,ksk)=sis;
        end
    end
end
xpop=Best_pop;
for i=1:1
    xin=xreg;
    n=len(i);
    if n==2
        neron(1,1:2)=Best_pop(i,1:2); 
    else
        layer=log(n)/log(2);
        if n==32
            layer=layer-1;
        end
        lay=1;
        neron(lay,:)=Best_pop(i,:);
        while lay<=(layer-1)
            qw=1;
            for oiu=1:2:len(i)-1
                neron(lay+1,qw)=neron(lay,oiu)*((Xin_col+1)^(2^(lay- 1)))+neron(lay,oiu+1);
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
            else
                for w=1:N_trn
                    A(w,:)=[1,xin(w,u),xin(w,v),xin(w,u)^2,xin(w,v)^2,xin(w,u)*xin(w,v)]; 
                end
                if metofcal==1
                    a=pinv(A'*A)*A'*(y(1:N_trn))';
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
                    a=CCc*BBb*AAa'*(y(1:N_trn))';
                    a=a';
                elseif metofcal==3
                    a=pinv(A'*A)*A'*(y(1:N_trn))';
                    a=a';
                    [x]=run_lm(xin,qwe,A,u,v,a)
                    a=x;
                elseif metofcal==4
                    a=[IG]';
                    [x]=run_lm(xin,qwe,A,u,v,a);
                    a=x; 
                end
                mutip(n_neroun,:)=a;
                n_neroun=n_neroun+1;
                clear CCc BBb AAa
                for w=1:cx
                    xnew(w,jj)=(a(1,1)+a(1,2)*xin(w,u)+a(1,3)*xin(w,v)+...
                        +a(1,4)*xin(w,u)^2+a(1,5)*xin(w,v)^2+a(1,6)*xin(w,u)*xin(w,v));
                end
            end
            jj=jj+1; 
        end
        len(i)=len(i)/2;
        clear xin
        xin=xnew;
        yout(1,:)=xnew(:,1)';
        clear xnew
    end
    if n==32
        for w=1:N_trn
            A(w,:)=[1,xin(w,1),xin(w,2),xin(w,1)^2,xin(w,2)^2,xin(w,1)*xin(w,2)]; 
        end
        if metofcal==1
            a=pinv(A'*A)*A'*(y(1:N_trn))';
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
            a=CCc*BBb*AAa'*(y(1:N_trn))';
            a=a';
        elseif metofcal==3
            a=pinv(A'*A)*A'*(y(1:N_trn))';
            a=a';
            [x]=run_lm(xin,qwe,A,u,v,a)
            a=x;
        elseif metofcal==4
              a=[IG]';
              [x]=run_lm(xin,qwe,A,u,v,a);
              a=x; 
        end
        clear CCc BBb AAa
        for w=1:cx
            xnew(w,1)=(a(1,1)+a(1,2)*xin(w,1)+a(1,3)*xin(w,2)+...
                +a(1,4)*xin(w,1)^2+a(1,5)*xin(w,2)^2+a(1,6)*xin(w,1)*xin(w,2)); 
        end
        xin=xnew

    end
    clear neron a xnew
     err(i)=0;
    err_predict=0;
    for we=1:qwe
        err(i)=err(i)+(xin(we)-y(we))^2; 
    end
    for we=qwe+1:qwe1
        err_predict(i)=err_predict(i)+(xin(we)-y(we))^2; 
    end
    yout=xin;
end
err_predict=err_predict/(qwe1-qwe);
err_Train=err/qwe;
mutip

