function [A] = UpdateLincomb_GD_emdbht(A,K,S,para,Obs)
% function solve min C1\sum_m \|Km-AmKmAm\|+ C2\|Am -\sum_lSlAl\|^2 + C3 \|A\|2_1
% A(i,:,m) : linear component for i th point in m th view 
%
% Input:
%

M=size(A,3);

for m=1:1:M
    [A,rtime(m),iteration(m)] = GD_emdbht(A,K,S,m,para,Obs);
    
end
end



function [Anew,spT,counter]=GD_emdbht(A,K,S,m,para,Obs)
alpha=1E-10;
beta=0.2;
eps=1E-6;
sigma = alpha;
Acur=A;

fcur=femdbht(Acur,K,S,para,Obs);
counter=1;


diff=10;

tstart=tic;

check=0;
while(diff>eps )

 
fold=fcur;
Aold=Acur;
 g=gemdbht(A,K,S,m,para,Obs);
 descent=g;
 [Acur,fcur]=getStepSize2(Aold,descent,g,K,S,m,para,sigma,beta,Obs);
 counter=counter+1;
 diff= (fold.T -fcur.T)/norm(fold.T+fcur.T);
 if (counter >10000  )
   Anew=Aold;
   spT= toc(tstart);

return;
  end


end
Anew=Aold;
   spT= toc(tstart);   
 
end





function [Anew,fnew]=getStepSize2(Aold,descent,g,K,S,m,para,sigma,beta,Obs)

step=1000; 
FLAG=1;
count=0;

fcur=femdbht(Aold,K,S,para,Obs);

gs=sigma*norm(g'*descent);

while(FLAG)
      count=count+1;
      temp=ProxL21(squeeze(Aold(:,:,m)-step*descent),para,step,size(Aold,3));
      Anew=Aold;
      Anew(:,Obs(m).id,m)=temp(:,Obs(m).id);
      unobs=setdiff([1:1:size(Aold,1)],Obs(m).id);
      Anew(:,unobs,m)=0;
      fnew =femdbht(Anew,K,S,para,Obs)  ;
    if(  fnew.T > fcur.T - 0.01*step*gs)
     step=beta*step;
     beta=beta*0.9;
    else
     return;
    end

    if norm(step*descent)<1E-18
            Anew=Aold;
            fnew=fcur;
            step=0;
       return;
    end
  end
end


           
  
  
