%%导入原始数据（光谱矩阵&理化值矩阵）
tic
load('RedBeautycitrus_650-1050_272_20211107.mat')
ROW_raw=size(X,1);
ValRatio=0.2;
TrNum=round(ROW_raw*(1-ValRatio));
data=ks_modification(X,Y,TrNum);
X=data.spec_train;
ROW_X=size(X,1);
Y=data.target_train;
X_val=data.spec_test;
Y_val=data.target_test;
%% 蒙特卡洛采样法剔除异常样本
%[X,Y]=OutlierElimination_MCS(X,Y,100);
%% 数据前处理-按照理化值Y降序排列
SD=std(Y, 1);
[row,col]=size(X);
Data=[Y X];
Data=sortrows(Data,-1);
Y=Data(:,1);
X=Data(:,2:(col+1));
X_raw=X;%原始光谱
FinalResult=[];
FinalEI=[];
num=0;
%% 波段选取-滑动窗口PLS部分参数设置
% min_interval=500;
% step_interval=100;
% window_interval=250;%50个变量间隔约为10个nm间隔
%测试用参数
% min_interval=500;
% step_interval=1000;
% window_interval=500;%50个变量间隔约为10个nm间隔
%% 样本集划分比例
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
                    SDc=std(Data_Label.target_train, 1);
                    SDp=std(Data_Label.target_test, 1);
                    SDv=std(Y_val);
                    CV = plscv(X,Y,20,10);
                    PLS=pls(Data_Label.spec_train,Data_Label.target_train,CV.optLV);
                    [ypred,Rp2,RMSEP]=plsval(PLS,Data_Label.spec_test,Data_Label.target_test,CV.optLV);
                    [yval_pred,Rv2,RMSEV]=plsval(PLS,X_val,Y_val,CV.optLV);
                    MPE=sum(abs(Y_val-yval_pred))/size(yval_pred,1);
                    MaxPE=max(abs(Y_val-yval_pred));
                    RPDp=SDp/RMSEP;
                    RPDv=SDv/RMSEV;
                    EI=2^(Rv2)*(0.5/RMSEV)*(2^(-(RMSEV/PLS.RMSEC-1)^2));
                    ModelingResult=[i,j,k,PLS.Rc2,PLS.RMSEC,Rp2,RMSEP,Rv2,RMSEV,MPE,MaxPE,RPDv,EI];
                    num=num+1;
                    FinalResult(num,:)=ModelingResult;
  FinalResult=sortrows(FinalResult,13,'descend'); 
  A=FinalResult;
  toc