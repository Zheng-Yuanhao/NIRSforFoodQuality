clc
clear

%%导入原始数据（光谱矩阵&理化值矩阵）
tic
load('citrus_2224_600_650-1050_channel1_20221207.mat')
ROW=size(X,1);
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
% min_interval=450;%MAYA光谱仪设定参数
% step_interval=90;
% window_interval=225;%45个变量间隔约为10个nm间隔
% min_interval=300;%辰昶光谱仪设定参数
% step_interval=60;
% window_interval=150;%30个变量间隔约为10个nm间隔
%测试用参数
min_interval=100;
step_interval=50;
window_interval=50;%50个变量间隔约为10个nm间隔
%% 样本集划分比例
for interval=col%min_interval:step_interval:col%%定义滑动窗口的范围
    for h=1%:window_interval:col-interval+1%%定义建模窗口的起始位置
        X=X_raw;
        X=X(:,h:h+interval-1);
        X_selection_raw=X;%记录波段选择后的原始光谱
        for i=1%1:3 %待使用的样本集划分比例情况种类
            if i==1
                TestRatio=0.2; %校正集：预测集=4:1；
            elseif i==2
                TestRatio=0.25; %校正集：预测集=3:1；
            elseif i==3
                TestRatio=1/3; %校正集：预测集=2:1；
            end
            %%样本预处理
            for j=6%1:6   %待使用的预处理方法数量
                    X=X_selection_raw;
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
                   %%样本划分
                    for k=1%:4 %待使用的样本划分方法数量
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
                    CV = plscv(Data_Label.spec_train,Data_Label.target_train,20,10);
                    PLS=pls(Data_Label.spec_train,Data_Label.target_train,CV.optLV);
                    [ypred,Rp2,RMSEP]=plsval(PLS,Data_Label.spec_test,Data_Label.target_test,CV.optLV);
                    waverange_a=Z(h);
                    waverange_b=Z(h+interval-1);
                    RPD=SDp/RMSEP;
                    CPI=3^(Rp2+(RPD-3)/3-abs(RMSEP/PLS.RMSEC-1));
                    ModelingResult=[waverange_a,waverange_b,h,interval,i,j,k,PLS.Rc2,PLS.RMSEC,Rp2,RMSEP,RPD,CPI];
                    num=num+1;
                    FinalResult(num,:)=ModelingResult;
                    end
            end
        end  
    end
end
  FinalResult=sortrows(FinalResult,13,'descend');
  A=FinalResult;
  toc