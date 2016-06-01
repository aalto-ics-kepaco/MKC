function [PredK,S,objective,iOutput]=MKCsdp(K,MID,para,init)

% this solve MKCsdp formulation 
% 
% input :
% K: kernels matrix size N x N x M. K(n1,n2,m) =kenrel values among n1 and n2 datapoints in m view
% MID: a matrix of NxM {0,1} matrix where MID(n,m) = 1 indicates nth data is known in mth view
% para: contain usedefined parameters 
%
%init : if init =1 then S is initialized by assignning all off-diagonal element with same values
%       otherwise it S is randomly initalized  
%
%Output:
% PredK =\hat(K) : predicted kernel matrices of size NxNxM
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
c1=para.c1;
c2=para.c2;
M=size(K,3);
N=size(K,1);

%idall=[1:1:N];

    diff= 10;
eps=1E-6;

for m=1:1:M
        Obs(m).id=find(MID(:,m)==1);
%        idall=intersect(idall,Obs(m).id);
end
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
%for m=1:1:M
%Kv(m,:)=reshape(K(idall,idall,m),1,numel(K(idall,idall,m)));
%end
%S=RelevantViews(Kv);

cvx_begin sdp 
variable PredK(N,N,M)
variable Loss(M)
variable Relevance(M)

minimize(c1*sum(Loss.*Loss) +c2*(sum(Relevance.*Relevance)))
for m=1:1:M
Loss(m)>=norm(squeeze(PredK(Obs(m).id,Obs(m).id,m))- squeeze(K(Obs(m).id,Obs(m).id,m)),'fro')
end
for m=1:1:M
   Relevance(m)>=norm(squeeze(PredK(:,:,m))- squeeze( sum(repmat(permute(S(m,:),[3,1,2]),[N,N,1]).* PredK,3)),'fro')
end
for m=1:1:M
    PredK(:,:,m) == semidefinite(N);
end

cvx_end
fcur = cvx_optval;
%fcur= 
iteration=0;
while diff >eps
for m=1:1:M
Kv(m,:)=reshape(PredK(:,:,m),1,numel(PredK(:,:,m)));
end
S=RelevantViews(Kv);
fold=fcur

PredK_old=PredK;
S_old=S;
iteration=iteration+1

iOutput(iteration).obj=fold;
iOutput(iteration).S=S_old;
iOutput(iteration).PredK=PredK_old;

cvx_begin sdp quiet
variable PredK(N,N,M)
variable Loss(M)
variable Relevance(M)

minimize(c1*sum(Loss.*Loss) +c2*(sum(Relevance.*Relevance)))
for m=1:1:M
Loss(m)>=norm(squeeze(PredK(Obs(m).id,Obs(m).id,m))- squeeze(K(Obs(m).id,Obs(m).id,m)),'fro')
end
for m=1:1:M
   Relevance(m)>=norm(squeeze(PredK(:,:,m))- squeeze( sum(repmat(permute(S(m,:),[3,1,2]),[N,N,1]).* PredK,3)),'fro')
end
for m=1:1:M
    PredK(:,:,m) == semidefinite(N);
end

cvx_end
fcur = cvx_optval;
diff=(fold-fcur)/norm(abs(fcur)+abs(fold))

end
PredK=PredK_old;
S=S_old;
objective=fold;

disp(PredK)
end
