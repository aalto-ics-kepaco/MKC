function[obj]=femdbht(A,K,S,para,Obs)
% this return object function value for MKCemdb(ht) formulation function solve min C1\sum_m \|Km-AmKmAm\|+ C2\|Am -\sum_lSmlAl\|^2 + C3 \|A\|2_1
% Output:
% obj.T = total objective values = para.c1 obj.K + para.c2 obj.Relevant +para.c3 obj.L21
% obj.K = Loss_within
% obj.L21 = L21 regularization
% obj.Relevant = Loss_between   

M=size(K,3); 
N=size(K,2);

obj.T=0;
obj.K=0;
%obj.CCA=0;
obj.Relevant=0;
obj.L21=0;

for l =1:1:M
temp(l).T=zeros(size(A,1),size(A,2));
for l2=1:1:M
	if l2==l
	else
           temp(l).T=S(l,l2)*A(:,:,l2);
	end
end
temp(l).D=A(:,:,l)*squeeze(K(:,:,l));
temp(l).hatk=temp(l).D*A(:,:,l)';
temp(l).B=K(:,:,l)-temp(l).hatk;
temp(l).E=temp(l).B(:,Obs(l).id)*temp(l).D(Obs(l).id,:);
temp(l).C=A(:,:,l)-temp(l).T;
end

for m=1:1:M
    obj.K=obj.K+((normS(K(Obs(m).id,Obs(m).id,m)-temp(m).hatk(Obs(m).id,Obs(m).id))^2)/(M*length(Obs(m).id)^2));
    for c=1:1:size(A,2)
           obj.L21=obj.L21+normS(A(:,c,m))/M;
     end	
  obj.Relevant=obj.Relevant+normS(temp(m).C)^2/(M*N);
end
obj.T=obj.K*para.c1+obj.Relevant*para.c2+obj.L21*para.c3;
end
