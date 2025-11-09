%%导入原始数据（光谱矩阵&理化值矩阵）
load('citrus_650-1050_600_20221109.mat')
SD=std(Y, 1);
[row,col]=size(X);
Data=[Y X];
Data=sortrows(Data,-1);
Y=Data(:,1);
X=Data(:,2:(col+1));
FinalResult=zeros(72,11);
FinalEI=zeros(72,4);
num=0;
i=2; %待使用的样本集划分比例情况种类
    if i==1
        TestRatio=0.2; %校正集：预测集=4:1；
    elseif i==2
        TestRatio=0.25; %校正集：预测集=3:1；
    elseif i==3
        TestRatio=0.3; %校正集：预测集=7:3；
    elseif i==4
        TestRatio=1/3; %校正集：预测集=2:1；
%   elseif i==5
%       TestRatio=0.4; %校正集：预测集=3:2;
    end
    CV = plscv(X,Y,20,10);
 j=6;
     %待使用的预处理方法数量
                    if j==1   %预处理为sg_smooth; 
                        [x_preprocessing] = sg_smooth(X,9,2,0);
                        X=x_preprocessing;  
                    elseif j==2  %预处理为normalize; 
                        [x_preprocessing] = normalize(X);
                        X=x_preprocessing;             
                    elseif j==3  %预处理为msc; 
                        [x_preprocessing] = msc(X);
                        X=x_preprocessing;
                    elseif j==4  %预处理为autoscales;        
                        [x_preprocessing] = autoscales(X);
                        X=x_preprocessing; 
                    elseif j==5  %预处理为meancenter;  
                        [x_preprocessing] = meancenter(X);
                        X=x_preprocessing;
                    elseif j==6  %预处理为snv;
                        [x_preprocessing] = snv(X);
                        X=x_preprocessing;
                    end
                k=2;%待使用的样本划分方法数量
                        TrainNum=round(row*(1-TestRatio));
                        if k==1   %样本划分方法为RS;
                            Data_Label=RandomSampleDivision(X,Y);
                        elseif k==2  %样本划分方法为KS;
                            Data_Label=ks_modification(X,Y,TrainNum);
                        elseif k==3  %样本划分方法为SPXY;
                            Data_Label=SPXY_modification(X,Y,TrainNum); 
                        elseif k==4  %样本划分方法为DOE;
                            Data_Label=descendingorder_extraction(X,Y,TestRatio);
                        end
                    %%PLS建模
                    SDp=std(Data_Label.target_test, 1);
                    PLS=pls(Data_Label.spec_train,Data_Label.target_train,CV.optLV);
                    [ypred,Rp2,RMSEP]=plsval(PLS,Data_Label.spec_test,Data_Label.target_test,CV.optLV);
                    Data_Label.target=[Data_Label.target_train;Data_Label.target_test];
                    [YPRED,R2,RMSE]=plsval(PLS,[Data_Label.spec_train;Data_Label.spec_test],[Data_Label.target_train;Data_Label.target_test],CV.optLV);
                    MPE=sum(abs(Data_Label.target_test-ypred))/size(ypred,1);
                    MaxPE=max(abs(Data_Label.target_test-ypred));
                    RPD=SDp/RMSEP;
                    EI=2^(Rp2)*(0.5/RMSEP)*(2^(-(RMSEP/PLS.RMSEC-1)^2));
                    ModelingResult=[i,j,k,PLS.Rc2,PLS.RMSEC,Rp2,RMSEP,MPE,MaxPE,RPD,EI];
                    num=num+1;
                    FinalResult(num,:)=ModelingResult;
                    A=sortrows(FinalResult,11,'descend');            