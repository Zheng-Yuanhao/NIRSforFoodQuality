%%导入原始数据（光谱矩阵&理化值矩阵）
tic
load('chucheng_600-1050_240_A1_20221108.mat')
ROW=size(X,1);
%%数据前处理-按照理化值Y降序排列
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
%%波段选取-滑动窗口PLS部分参数设置
%%样本集划分比例
for ii=1:col-1
    for jj=ii+149:col 
        X=X_raw;
        X=X(:,ii:jj);
        X_selection_raw=X;%记录波段选择后的原始光谱
        for i=1:3  %待使用的样本集划分比例情况种类
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
            for j=1:6   %待使用的预处理方法数量
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
                    for k=1:4 %待使用的样本划分方法数量
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
                    %CV = plscv(X,Y,20,10);
                    CV = plscv(Data_Label.spec_train,Data_Label.target_train,20,10);
                    PLS=pls(Data_Label.spec_train,Data_Label.target_train,CV.optLV);
                    [ypred,Rp2,RMSEP]=plsval(PLS,Data_Label.spec_test,Data_Label.target_test,CV.optLV);
                    MPE=sum(abs(Data_Label.target_test-ypred))/size(ypred,1);
                    MaxPE=max(abs(Data_Label.target_test-ypred));
                    Outlier_num=ROW-row;
                    waverange_a=Z(ii);
                    waverange_b=Z(jj);
                    RPD=SDp/RMSEP;
                    EI_0=2^(Rp2)*(0.5/RMSEP)*(2^(-abs((RMSEP/PLS.RMSEC-1))));
                    EI_1=2^(Rp2+(RPD-2)-abs(RMSEP/PLS.RMSEC-1));
                    EI_2=2^(Rp2+(RPD/2-1)-abs(RMSEP/PLS.RMSEC-1));
                    ModelingResult=[waverange_a,waverange_b,i,j,k,PLS.Rc2,PLS.RMSEC,Rp2,RMSEP,MPE,MaxPE,RPD,EI_0,EI_1,EI_2];
                    num=num+1;
                    FinalResult(num,:)=ModelingResult;
                    end
            end
        end  
    end
end
  FinalResult=sortrows(FinalResult,12,'descend');
  A=FinalResult;
  toc