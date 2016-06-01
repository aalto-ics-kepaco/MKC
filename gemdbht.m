function [grad] = gemdbht(A,K,S,m,para,Obs,Kc)
% this calculate gradient step of objective function for MKC_emdb(ht) with respect to A^{(m)}
% Output :
% grad :  (\partial f / \partial A^{(m)})

grad=zeros(size(A,1),size(A,2));
eps=1E-10;
M=size(K,3); 
N=size(K,2);
  
for l =1:1:M
temp(l).T=zeros(size(A,1),size(A,2));
for l2=1:1:M
	if l2==l
	else
           temp(l).T=S(l,l2)*A(:,:,l2);
	end
end

temp(l).D=A(:,:,l)*K(:,:,l);
temp(l).hatk=temp(l).D*A(:,:,l)';
temp(l).B=K(:,:,l)-temp(l).hatk;
temp(l).E=temp(l).B(:,Obs(l).id)*temp(l).D(Obs(l).id,:);
temp(l).C=A(:,:,l)-temp(l).T;
end

 Part1=zeros(N,N);


Part1=-temp(m).E;
grad= grad +4*para.c1*Part1/(M*length(Obs(m).id)^2);
                          
part= zeros(size(A,1),size(A,2));
for l =1:1:M
if l==m
   part=part+temp(m).C;
else
part=part-S(l,m)*temp(l).C;
end
end
grad= grad + 2*para.c2*part/(M*N);


clear temp part Part1;
             
    return;     
end














































