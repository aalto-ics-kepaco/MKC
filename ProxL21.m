function [delA] = ProxL21(delA,para,stepsize,M)
% this calculate proximal operator on L21 norm
% Output :
% grad :  Prox(grad)

for c=1:1:size(delA,2)
	if normS(delA(:,c))>eps
		delA(:,c)=max(0, (1- stepsize*para.c3/(M*normS(delA(:,c))))*delA(:,c));
	end
end

    return;     
end














































