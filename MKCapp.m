function [PredK,A,S,obj,iOutput] = MKCapp(K,MID,para,init)

% this solve MKCkernel formulation 
% 
% input :
% K: kernels matrix size N x N x M. K(n1,n2,m) =kenrel values among n1 and n2 datapoints in m view
% MID: a matrix of NxM {0,1} matrix where MID(n,m) = 1 indicates nth data is known in mth view
% para: contain usedefined parameters 
%
%init : if init =1 then A^{(m)} is initialized by non diagonal selfrepresentative matrix of K^{(m)}
%       otherwise it A^{(m)} is randomly initalized  
%
%Output:
% PredK =\hat(K) : predicted kernel matrices of size NxNxM
% A : learnt reconstruction matrix
% S : Inter-view similarity matrix
% objective values
% iOutput : intermediate objective function 

% (c) Sahely Bhadra
% sahely.bhadra@aalto.fi
%
% Algorithm is explained as in the linked document
% http://arxiv.org/abs/1602.02518
%
% Jun. 1, 2016.
 
  

tstart=tic;
M=size(K,3);
N=size(K,2);

diff= 10;
eps=1E-5;

for m=1:1:M
        Obs(m).id=find(MID(:,m)==1);
end

% initialization 
if init==1
S=ones(M,M);
S=S-diag(diag(S));
    for m=1:1:M
        S(m,:)=ProjSimplex(S(m,:));
    end
    S=S-diag(diag(S));
else
    S=randn(M,M);
    S=S-diag(diag(S));
    for m=1:1:M
        S(m,:)=ProjSimplex(S(m,:));
    end
    S=S-diag(diag(S));
end

for m=1:1:M 
   if init==1
         A(:,:,m)=zeros(size(K(:,:,m)));
         A(:,Obs(m).id,m)=randn(N,length(Obs(m).id));
         A(Obs(m).id,Obs(m).id,m)=RelevantViews(K(Obs(m).id,Obs(m).id,m));
   else
       A(:,:,m)=zeros(size(K(:,:,m)));
       A(:,Obs(m).id,m)=randn(N,length(Obs(m).id));
   end
end
   
fcur=fkernel(A,K,S,para,Obs);
iteration=0;

% parametes are optimized by block wise decsend 

while(diff>eps)
 
 Sold=S;
 Aold=A;
 fold=fcur;
iteration=iteration+1;
iOutput.intermediate(iteration).f=fold; % stoer objecttive function value in each iteration 
iOutput.iteration=iteration;

    %updating A
    A = UpdateLincomb_GD_kernel(Aold,K,Sold,para,Obs);
   
    %for m=1:1:M 
    %L(m,:)=reshape(A(:,:,m),1,numel(A(:,:,m)));
    %end
    for m=1:1:M
    tempK = A(:,:,m)*K(:,:,m)*A(:,:,m)'; 
    L(m,:)=reshape(tempK,1,numel(tempK));
    end
    S=RelevantViews(L);
    fcur=fkernel(A,K,S,para,Obs);
    iteration =iteration +1;
    diff=(fold.T-fcur.T)/norm(abs(fcur.T)+abs(fold.T))
     if (iteration >10000  )
            A=Aold;
           
            S=Sold;
            obj=fkernel(A,K,S,para,Obs);
            spT= toc(tstart);
              for m=1:1:M
            PredK(:,:,m)=squeeze(A(:,:,m))*squeeze(K(:,:,m))*squeeze(A(:,:,m))';
            end
    return;
     end
  
end
            A=Aold;
            S=Sold;
            obj=fkernel(A,K,S,para,Obs);
            spT= toc(tstart);
            for m=1:1:M
            PredK(:,:,m)=squeeze(A(:,:,m))*squeeze(K(:,:,m))*squeeze(A(:,:,m))';
            end
end
