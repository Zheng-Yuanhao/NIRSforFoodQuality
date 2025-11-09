function [R2,RMSE,RPD]=R2RMSE_calculation(ytest,ypred)
true_y=ytest;
predict_y=ypred;
tempdata=(true_y-predict_y).^2;
tempdata2=(true_y-mean(true_y)).^2;
R2=1-(sum(tempdata)/sum(tempdata2));
RMSE=sqrt(sumsqr(ypred-ytest)/length(ytest));
RPD=std(ytest)/RMSE;
end