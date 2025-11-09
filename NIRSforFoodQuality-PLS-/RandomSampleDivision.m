function Data_Label=RandomSampleDivision(X,Y)
data=X;
label=Y;
% 测试数据占全部数据的比例
testRatio = 0.2;
% 训练集索引
trainIndices = crossvalind('HoldOut', size(data, 1), testRatio);
% 测试集索引
testIndices = ~trainIndices;
% 训练集和训练标签
spec_train = data(trainIndices, :);
target_train = label(trainIndices, :);
% 测试集和测试标签
spec_test = data(testIndices, :);
target_test = label(testIndices, :);
Data_Label.spec_train=spec_train;
Data_Label.target_train=target_train;
Data_Label.spec_test=spec_test;
Data_Label.target_test=target_test;
end