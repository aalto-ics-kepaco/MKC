function [Kn] = normK(K)

for m=1:1:size(K,3)
Kn(:,:,m)=K(:,:,m)./(sqrt(diag(K(:,:,m))*diag(K(:,:,m))'));
end

end
