function [ypred,Rp2,RMSEP]=plsval(plsmodel,Xtest,ytest,nLV)
%+++ Compute prediction errors on a test set
%+++ plsmodel: a structural data obtained from the function pls.m in this directory;
%+++ Xtest: test samples
%+++ ytest: y values for test samples (for "really" new samples, no y values available
%+++ nLV: number of latent variables for calibration models.
%+++ Hongdong Li,Oct.21,2007;
if nargin<4;nLV=size(plsmodel.X_scores,2);end;
Xtest=[Xtest ones(size(Xtest,1),1)]; 
ypred=Xtest*plsmodel.regcoef_original_all(:,nLV);
true_y=ytest;
predict_y=ypred;
tempdata=(true_y-predict_y).^2;
tempdata2=(true_y-mean(true_y)).^2;
Rp2=1-(sum(tempdata)/sum(tempdata2));
RMSEP=sqrt(sumsqr(ypred-ytest)/(length(ytest)-1));






