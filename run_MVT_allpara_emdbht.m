
function run_MVT_allpara_emdbht(c,name)
%load(['../Toyresult/rmse_MFEAT_small2.mat'],'rmse','para');


%count=size(para,1)
%str=count +1;
count=0;
str=count+1;
%for num=[1,10]
for i=[1000]
    c1=i;
    for j=[1, 10]
        c2=j; %relevence
        for k=[ 0.0001 0.001 0.01 0.1 1 10 ]
            c3=k; %l21
            %for c4=[1 10 100]  %cca
            			count=count+1;
            			para(count,:)=[c1,c2,c3];
            %        end
        end
     end
end
%end
stp=count
count=c
c1=para(count,1);
c2=para(count,2);
c3=para(count,3);
inpara=[c1,c2,c3];
 [rmse]=MKCemdbht_run(name,3,c1,c2,c3);
fname =[name,'_emdbht'];
system(['mkdir ../',fname]);
save(['../',fname,'/rmse',num2str(c),'.mat'],'rmse','para','c','inpara');
   


end

