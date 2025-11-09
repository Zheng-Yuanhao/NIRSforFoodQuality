%%此程序为matlab编程实现的BP神经网络
% 清空环境变量
tic
clear
close all
clc
num=5;
Distance=cell(1,num);
Ans=cell(1,num);
%%第一步 读取数据
%load('mango_684-1990_10243_2015-2019-bp.mat')
load('pharmaceuticaltable-602-600-1622-2002-bp.mat')
%load('corn_oil_1100-2498_80_m5-bp.mat')
%%第二步 设置训练数据和预测数据
input_train = input_train';
output_train =output_train';
input_test = input_test';
output_test =output_test';
%节点个数
inputnum=size(input_train,1); % 输入层节点数量n
outputnum=size(output_train,1); % 输出层节点数量m
a=10;%调节常数1-10
hiddennum=ceil(sqrt(inputnum+outputnum)+a);% 隐含层节点数量一般位于sqrt(m+n)+a，或log2(n)取整，必须小于输入层节点n-1,
%%第三本 训练样本数据归一化
[inputn,inputps]=mapminmax(input_train);%归一化到[-1,1]之间，inputps用来作下一次同样的归一化
[outputn,outputps]=mapminmax(output_train);
for i=1:num
        %%第四步 构建BP神经网络
        net=newff(inputn,outputn,hiddennum,{'tansig','purelin'},'trainlm');% 建立模型，传递函数使用purelin，采用梯度下降法训练
        net.trainParam.showWindow = false;
        net.trainParam.showCommandLine = false; 
        W1= net. iw{1, 1};%输入层到中间层的权值
        B1 = net.b{1};%中间各层神经元阈值
        W2 = net.lw{2,1};%中间层到输出层的权值
        B2 = net. b{2};%输出层各神经元阈值
        %%第五步 网络参数配置（训练次数，学习速率，训练目标最小误差等
        anss=[];
        for j=1:100
        net.trainParam.epochs=1;         % 训练次数，这里设置为100次
        net.trainParam.lr=0.01;                   % 学习速率，这里设置为0.01
        net.trainParam.goal=0.00001;                    % 训练目标最小误差，这里设置为0.00001
        %%第六步 BP神经网络训练
        net=train(net,inputn,outputn);%开始训练，其中inputn,outputn分别为输入输出样本
        %%第七步 测试样本归一化
        inputn_test=mapminmax('apply',input_test,inputps);% 对样本数据进行归一化
        %%第八步 BP神经网络预测
        an_train=sim(net,inputn); %用训练好的模型对训练集进行仿真
        an_test=sim(net,inputn_test); %用训练好的模型对测试集进行仿真
        %%第九步 预测结果反归一化与误差计算  
        train_simu=mapminmax('reverse',an_train,outputps); %把训练集仿真得到的数据还原为原始的数量级
        train_error=train_simu-output_train;      %训练集预测值和真实值的误差
        [Rc2,RMSEC,RPDc]=R2RMSE_calculation(train_simu,output_train);
        test_simu=mapminmax('reverse',an_test,outputps); %把测试集仿真得到的数据还原为原始的数量级
        test_error=test_simu-output_test;      %测试集预测值和真实值的误差
        [Rp2,RMSEP,RPDp]=R2RMSE_calculation(test_simu,output_test);
        EI_0=3^(Rp2+(RPDp-3)/3-(abs(RMSEP/RMSEC-1)));
        anss(j,:)=[Rp2 (RPDp-3)/3 EI_0];
                if j==1
                   distance=0;
                else
                   distance=distance+sqrt((anss(j,1)-anss(j-1,1))^2+(anss(j,2)-anss(j-1,2))^2);    
                end    
        end
Distance{i}=[Distance{i},distance];
Ans{i}=[Ans{i},anss];
end
toc
%%第十步 真实值与预测值误差比较
% figure('units','normalized','position',[0.119 0.2 0.38 0.5])
% plot(output_test,'bo-')
% hold on
% plot(test_simu,'r*-')
% hold on
% plot(test_error,'square','MarkerFaceColor','b')
% legend('期望值','预测值','误差')
% xlabel('数据组数')
% ylabel('样本值')
% title('BP神经网络测试集的预测值与实际值对比图')
% [c,l]=size(output_test);
% MAE1=sum(abs(test_error))/l;
% MSE1=test_error*test_error'/l;
% RMSE1=MSE1^(1/2);

% disp(['-----------------------误差计算--------------------------'])
% disp(['隐含层节点数为',num2str(hiddennum),'时的误差结果如下：'])
% disp(['平均绝对误差MAE为：',num2str(MAE1)])
% disp(['均方误差MSE为：       ',num2str(MSE1)])
% disp(['均方根误差RMSE为：  ',num2str(RMSE1)])
% 附
% eval(['web ', char([104	116	116	112	115	58	47	47	98	108	111	103	46	99	115	100	110	46	110	101	116	47	113	113	95	53	55	57	55	49	52	55	49	47	97	114	116	105	99	108	101	47	100	101	116	97	105	108	115	47	49	50	49	55	54	55	48	48	52	32	45	98	114	111	119	115	101	114])])
% eval(['web ', char([104,116,116,112,115,58,47,47,109,98,100,46,112,117,98,47,111,47,98,114,101,97,100,47,109,98,100,45,89,90,109,84,109,112,116,118,32,45,98,114,111,119,115,101,114])])
% eval(['web ', char([104,116,116,112,115,58,47,47,109,98,100,46,112,117,98,47,111,47,117,112,115,95,100,111,119,110,115,47,119,111,114,107,32,45,98,114,111,119,115,101,114])]) 