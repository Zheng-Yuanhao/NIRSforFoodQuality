%%导入原始数据（光谱矩阵&理化值矩阵）
tic
load('RedBeautycitrus_650-1050_450_20221103.mat')
load('Redcitrus_650-1050_390_20221110_val1.mat')
%% 蒙特卡洛采样法剔除异常样本
% [X,Y]=OutlierElimination_MCS(X,Y,100);
%% 数据前处理-按照理化值Y降序排列
SD=std(Y, 1);
[row,col]=size(X);
Data=[Y X];
Data=sortrows(Data,-1);
Y=Data(:,1);
X=Data(:,2:(col+1));
X_raw=X;
FinalResult=zeros(72,12);
FinalEI=zeros(72,5);
num=0;
%% 样本集划分比例
for i=1 %待使用的样本集划分比例情况种类
    if i==1
        TestRatio=0.2; %校正集：预测集=4:1；
    elseif i==2
        TestRatio=0.25; %校正集：预测集=3:1；
%     elseif i==3
%         TestRatio=0.3; %校正集：预测集=7:3；
    elseif i==3
        TestRatio=1/3; %校正集：预测集=2:1；
%   elseif i==5
%       TestRatio=0.4; %校正集：预测集=3:2;
    end
    %%样本预处理
    for j=6%待使用的预处理方法数量
            X=X_raw;
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
            for k=3 %待使用的样本划分方法数量
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
            MPE=sum(abs(Data_Label.target_test-ypred))/size(ypred,1);
            MaxPE=max(abs(Data_Label.target_test-ypred));
            RPD=SDp/RMSEP;
            EI=2^(Rp2)*(0.5/RMSEP)*(2^(-(RMSEP/PLS.RMSEC-1)^2));
            ModelingResult=[i,j,k,PLS.Rc2,PLS.RMSEC,Rp2,RMSEP,MPE,MaxPE,CV.optLV,RPD,EI];
            num=num+1;
            FinalResult(num,:)=ModelingResult;
            end
    end
end
  FinalResult=sortrows(FinalResult,12,'descend'); 
  A=FinalResult;
  %%模型参数检验
  a_name='VERSNVPLS';
  a_number=col;
  a_bias= mean(Data_Label.target_train);
  a_mean = mean(Data_Label.spec_train,1);
  a_mean = a_mean';
  Mode = [(1:a_number)'-1+a_number,a_mean, PLS.regcoef_original(1:end-1,:)];
  %Mode = [(1:a_number)'-1+num_wave,a_mean, PLS.regcoef_original(1:end-1,:)];
  mm=sum((Data_Label.spec_train-[Mode(:,2)]').*[Mode(:,3)]',2)+a_bias;
  %mm应该等于PLS.y_fit
  %%外部数据模型验证
  X_val=snv(X_val);
  [yval_pred,Rv2,RMSEV]=plsval(PLS,X_val,Y_val,CV.optLV);
  toc