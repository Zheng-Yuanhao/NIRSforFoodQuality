load Channel1Data.mat
X_raw=X;
Result=zeros(16,7);
%波段1：650-900nm；波段2：650-910nm；波段3：650-920nm...，波段16：650-1050nm
for i=1:16
    if i<3 
        m=1149+48*(i-1);%650-910nm
    elseif i<9 
        m=1245+49*(i-3);%650-970nm
    elseif i<13
        m=1540+50*(i-9);%650-1020nm
    else
        m=1741+51*(i-13);%650-1050nm          
    end      
    X=X(:,1:m);
    CV = plscv(X,Y,20,10,'center');
    Data_Label=RandomSampleDivision(X,Y);
    %Data_Label=SPXY_modification(X,Y,88);
    PLS=pls(Data_Label.spec_train,Data_Label.target_train,CV.optLV);
    % load Data_all.mat
    % PLS=pls(xtrain,ytrain,11);
    [ypred,Rp2,RMSEP]=plsval(PLS,Data_Label.spec_test,Data_Label.target_test,CV.optLV);
    %[ypred,Rp2,RMSEP]=plsval(PLS,Xtest,ytest,11);
%     scatter(Data_Label.target_test,ypred);
%     hold on
%     xx=7:0.05:15;
%     y=xx;
%     plot(xx,y,'k','LineWidth',2)
    MPE=sum(abs(Data_Label.target_test-ypred))/size(ypred,1);
    MaxPE=max(abs(Data_Label.target_test-ypred));
    ModelingResult=[PLS.Rc2,PLS.RMSEC,Rp2,RMSEP,MPE,MaxPE,CV.optLV];
    Result(i,:)=ModelingResult;
    X=X_raw;
end