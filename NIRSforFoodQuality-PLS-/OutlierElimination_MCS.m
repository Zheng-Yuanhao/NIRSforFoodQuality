function [X,Y]=OutlierElimination_MCS(X,Y,N)
%+++ Monte-Carlo Sampling method for Outlier Elimination.
%+++ X: sample matrix: M x P.
%+++ Y: Response variable: p-dimensional.
%+++ N: Number of Monte Carlo Sampling
tic
A=15;
ratio=0.8;
OPT=1;
[Mx,Nx]=size(X);
A=min([size(X) A]);
nc=floor(Mx*ratio);
nv=Mx-nc;
yp=[];yc=[];
error_pre=[];
y=[];
%%%确定超参数个数
for i=1:N
    % generate sub-training set and sub-test set
    index=randperm(Mx);
    calk=index(1:nc);
    testk=index(nc+1:Mx);
    Xcal=X(calk,:);ycal=Y(calk);
    Xtest=X(testk,:);ytest=Y(testk);
    % data pretreatment
    Xcal=snv(Xcal);
    Xtest=snv(Xtest);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLS modeling 
    CV=plscv(Xcal,ycal,20,10);
    myPLS_vis=pls(Xcal,ycal,CV.optLV,'center');
    [ytest_p,RMSEP_test,R2_test]=plsval(myPLS_vis,Xtest,ytest,CV.optLV);
    error_pre(i,testk)=ytest_p-ytest; 
    if OPT==1;if mod(i,100)==0;fprintf('The %dth sampling for MCCV finished.\n',i);end;end
end
% calculate the mean and std of predictive error for N 
for i=1:size(error_pre,2)
    Sub=[];
    Mean(i)=sum(error_pre(:,i))/sum(error_pre(:,i)~=0);
    ind=find(error_pre(:,i));
    Sub=error_pre(ind,i)-Mean(i);
    Std(i)=sum(Sub.^2)/(sum(error_pre(:,i)~=0)-1);
end
% calculate the mean and std of predictive error for each sample
thre_mean_Mean=mean(Mean);
thre_Std_Mean=std(Mean);
thre_mean_Std=mean(Std);
thre_Std_Std=std(Std);
% 3 sigma
Outliers_index_Mean=find(Mean>=(thre_mean_Mean+3*thre_Std_Mean) | Mean<=(thre_mean_Mean-3*thre_Std_Mean));
Outliers_index_Std=find(Std>=(thre_mean_Std+3*thre_Std_Std) | Std<=(thre_mean_Std-3*thre_Std_Std));
Outliers_index=[Outliers_index_Mean,Outliers_index_Std];UnOutliers_index=unique(Outliers_index);
figure(1)
scatter(Mean,Std)
hold on
yy=0:max(Std)/100:max(Std);
xx=(thre_mean_Mean-3*thre_Std_Mean)*ones(1,101);
x=(thre_mean_Mean+3*thre_Std_Mean)*ones(1,101);
xxx=min(Mean):((max(Mean)-min(Mean))/100):max(Mean);
yyy=(thre_mean_Std-3*thre_Std_Std)*ones(1,101);
yyyy=(thre_mean_Std+3*thre_Std_Std)*ones(1,101);
plot(xx,yy)
hold on
plot(x,yy)
hold on
plot(xxx,yyy)
hold on
plot(xxx,yyyy)
% figure(1)
% plot(X(Outliers_index_Mean,:)')
% figure(2)
% plot(X(Outliers_index_Std,:)')
% figure(3)
% plot(X(UnOutliers_index,:)')
X(UnOutliers_index,:)=[];Y(UnOutliers_index)=[];
toc;