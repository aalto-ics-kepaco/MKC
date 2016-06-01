function [rmse]=MKCemdnht_run(name,Rep,c1,c2,c3)

para.c1=c1;
para.c2=c2;
para.c3=c3;
fname=[name,'_emdbht'];
system(['mkdir ../',fname]);
str=[num2str(c1),'_',num2str(c2),'_',num2str(c3)];

load(['../data/',name,'.mat'])
M=size(K,3);
if length(size(IDmisTest))==3
MISS=size(IDmisTest,3);
end
if length(size(IDmisTest))==4
MISS=size(IDmisTest,4);
end

for nInd=1:1:MISS

for r=1:1:Rep
     for cv=2:1:2
         if cv ==1   
            idts = squeeze([IDtr(r,:),IDts(r,:)]);
            Kts=K(idts,idts,:);
            MID=ones(length(IDtr(r,:)),M);
	if length(size(IDmisTest))==3
                MID=[MID; squeeze(IDmisTest(:,:,nInd))];
            end
            if length(size(IDmisTest))==4
                MID=[MID; squeeze(IDmisTest(r,:,:,nInd))];
            end
         else
            idts = squeeze([IDtr(r,:),IDval(r,:)]);
            Kts=K(idts,idts,:);
            MID=ones(length(IDtr(r,:)),M);
            if length(size(IDmisCv))==3
                MID=[MID; squeeze(IDmisCv(:,:,nInd))];
            end
            if length(size(IDmisCv))==4
                MID=[MID; squeeze(IDmisCv(r,:,:,nInd))];
            end   

        end


      if exist(['../',fname,'/sF',name,'miss_',num2str(nInd),'r_',num2str(r),'_cv_',num2str(cv),'_',str,'Model.mat'])==2

        disp(['../',fname,'/sF',name,'miss_',num2str(nInd),'r_',num2str(r),'_cv_',num2str(cv),'_',str,'Model.mat'])
        load(['../',fname,'/sF',name,'miss_',num2str(nInd),'r_',num2str(r),'_cv_',num2str(cv),'_',str,'Model.mat'],'Mer');
      
      else
	        O=9999999999999;
                        
        	for init=1:1:50
                         stime = tic;        
              		 [tPredK,tA,tS,tobjective,tiOutput]=MKCemdbht(Kts,MID,para,init);
                         runtime(init)=toc(stime);
               		if tobjective.T <O
                   	   PredK=tPredK;
	                   Model.A=tA;
        	           Model.S=tS;
                	   Model.objective=tobjective;
                           Model.iOutput=tiOutput;
                           O=tobjective.T
               		end
                       
            	end

       
	        Model.Sall=Model.S;
      		Model.S(find(Model.S<=0.05))=0;
          	Model.Aall=Model.A
      		Model.A(find(Model.A<=0.001))=0;
      		M=size(K,3)
      		for m=1:1:M
          		unObs(m).id=find(MID(:,m)==0);
          		if length(unObs(m).id)>0
				info(m).Er=diag((PredK(unObs(m).id,:,m)-Kts(unObs(m).id,:,m))*(PredK(unObs(m).id,:,m)-Kts(unObs(m).id,:,m))')./diag(Kts(unObs(m).id,:,m)*Kts(unObs(m).id,:,m)');
          		else
				info(m).Er=0;
                        end
			Mer(m)=mean(info(m).Er);
          		Ser(m)=std(info(m).Er);
          		info(m).core=find(sum(abs(Model.A(:,:,m)))>0.3)
 			Ncore(m)=length(info(m).core);
      		end
                save(['../',fname,'/sF',name,'miss_',num2str(nInd),'r_',num2str(r),'_cv_',num2str(cv),'_',str,'Model.mat'] ,'Model','PredK','Kts','MID','unObs','info','Ncore','Mer','Ser','runtime');
	end
      	rmse(r,nInd,cv)=mean(Mer);
      end
      end
   end
return;
end


