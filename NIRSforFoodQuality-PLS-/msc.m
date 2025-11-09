%% 多元散射校正 Multivariate Scatter Correction
function [x_preprocessing] = msc(x)
% Multiplicative Scatter Correction
% [x_msc] = msc(x)
% input:
% x (samples x variables) data to preprocess
% output:
% x_msc (samples x variables) preprocessed data
 [m,n]=size(x);
 me=mean(x);
 for i=1:m
     p=polyfit(me,x(i,:),1);
     x_preprocessing(i,:)=(x(i,:)-p(2)*ones(1,n))./(p(1)*ones(1,n));
 end