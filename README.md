# MKC
Multi-view Kernel Completion

(c) Sahely Bhadra
sahely.bhadra@aalto.fi
Details of the software are available in http://arxiv.org/abs/1602.02518

Jun. 1, 2016.
This package contains  following three version of proposed multi-view kernel completion method along with supporting function and scripts

1. MKCsdp
2. MKCapp
3. MKCemdb(ht)


Details descriare :




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

function [PredK,A,S,obj,iOutput] = MKCemdbht(K,MID,para,init)

% this solve MKCembdht formulation 
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

