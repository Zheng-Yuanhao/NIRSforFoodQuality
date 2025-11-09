%% 归一化 Normalize
function [x_preprocessing] = normalize(x)
% Normalize matrix rows dividing by its norm
% [nx] = normaliz(x)
% input:
% x (samples x variables)   data to normalize
% output:
% nx (samples x variables)  normalized data
% By Cleiton A. Nunes
% UFLA,MG,Brazil
[m,n]=size(x);
x_preprocessing=x;
nm=zeros(m,1);
for i = 1:m
nm(i)=norm(x_preprocessing(i,:));
x_preprocessing(i,:)=x_preprocessing(i,:)/nm(i);
end