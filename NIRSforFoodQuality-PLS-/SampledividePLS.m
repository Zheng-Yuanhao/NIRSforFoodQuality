%% 导入原始数据（光谱矩阵&理化值矩阵）
load('Twochannelorange_channel1.mat')
SD=std(Y, 1);
[row,col]=size(X);
Data=[Y X];
Data=sortrows(Data,-1);
Y=Data(:,1);
X=Data(:,2:(col+1));
%% 交叉验证（无预处理）
CV = plscv(X,Y,20,10);
%% 预处理方法
[x_preprocessing] = meancenter(X);
X=x_preprocessing;
%% 样本划分
%Data_Label=RandomSampleDivision(X,Y);
Data_Label=ks_modification(X,Y,88);
%% PLS建模
PLS=pls(Data_Label.spec_train,Data_Label.target_train,CV.optLV);
[ypred,Rp2,RMSEP]=plsval(PLS,Data_Label.spec_test,Data_Label.target_test,CV.optLV);
MPE=sum(abs(Data_Label.target_test-ypred))/size(ypred,1);
MaxPE=max(abs(Data_Label.target_test-ypred));
ModelingResult=[PLS.Rc2,PLS.RMSEC,Rp2,RMSEP,MPE,MaxPE,CV.optLV];
%% 预测结果绘图
scatter(Data_Label.target_test,ypred);
hold on
xx=7:0.05:15;
y=xx;
plot(xx,y,'k','LineWidth',2)