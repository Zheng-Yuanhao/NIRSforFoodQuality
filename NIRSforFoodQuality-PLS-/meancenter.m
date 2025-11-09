%% 数据中心化(Mean centering)
function [x_preprocessing] = meancenter(x)
% Mean center scales matrix to zero mean
% [mcx,mx] = center(x)
% input:
% x 	data to mean center
% output:
% ax	mean center data
% mx	means of data
% By Cleiton A. Nunes
% UFLA,MG,Brazil
[m,n]=size(x);
mx=mean(x);
x_preprocessing=(x-mx(ones(m,1),:));