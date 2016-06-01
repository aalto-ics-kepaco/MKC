function[obj]=fkernel(A,K,S,para,Obs)
% this return object function value for MKCkernel formulation
% Output:
% obj.T = total objective values = para.c1 obj.K + para.c2 obj.Relevant +para.c3 obj.L21
% obj.K = Loss_within
% obj.L21 = L21 regularization
% obj.Relevant = Loss_between   

M=size(K,3); 
N=size(K,2);

obj.T=0;
obj.K=0;
obj.Relevant=0;
obj.L21=0;
for m=1:1:M
   temp(m).hatk= A(:,:,m)*K(:,:,m)*A(:,:,m)';
end

for m=1:1:M

    obj.K=obj.K+((normS(K(Obs(m).id,Obs(m).id,m)-temp(m).hatk(Obs(m).id,Obs(m).id))^2)/(M*length(Obs(m).id)^2));

                for c=1:1:size(A,2)
                	obj.L21=obj.L21+normS(A(:,c,m))/M;
           	    end	
    temp2=zeros(size(A,1),size(A,2));

     for l=1:1:size(K,3)
            if l==m
            else
                temp2 =temp2 + S(m,l)*temp(l).hatk; % relevence part
            end
     end
obj.Relevant=obj.Relevant+normS(temp(m).hatk-temp2)^2/(M*N);
    
end

obj.T=obj.K*para.c1+obj.Relevant*para.c2+obj.L21*para.c3;

end
