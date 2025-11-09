function stop=plotTrainingAccuracy(info)
    % 在每个epoch之后输出训练准确率和测试准确率
    if info.State == "end"
          %disp(fieldnames(info.TrainingManager.UserData))
          disp(1)
          %         %预测结果
%         train_YPred = predict(net,trainD);%训练集的预测值
%         test_YPred = predict(net,testD);%预测集的预测值
%         train_YPred=double(train_YPred');
%         test_YPred=double(test_YPred');%输出是n*1的single型数据，要转换为1*n的double是数据形式
%         %反归一化
%         train_predict_value=method('reverse',train_YPred,output_ps);train_predict_value=double(train_predict_value);
%         train_true_value=method('reverse',targetD_train,output_ps);train_true_value=double(train_true_value);
%         test_predict_value=method('reverse',test_YPred,output_ps);test_predict_value=double(test_predict_value);
%         test_true_value=method('reverse',targetD_test,output_ps);test_true_value=double(test_true_value);
%         %%模型评价
%         [Rc2,RMSEC,RPDc]=R2RMSE_calculation(train_predict_value,train_true_value);
%         [Rp2,RMSEP,RPDp]=R2RMSE_calculation(test_predict_value,test_true_value);
%         EI_0=3^(Rp2+(RPDp-3)/3-(abs(RMSEP/RMSEC-1)));
%         anss=[Rp2 (RPDp-3)/3 EI_0];
%         info.anss=anss;
%         % 输出训练结果和测试结果
%         fprintf('Epoch %d: Training RMSE=%.4f',info.Epoch,EI_0)
    end
end