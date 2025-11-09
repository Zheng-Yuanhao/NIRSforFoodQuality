function Data_Label=descendingorder_extraction(X,Y,TestRatio)
row=size(X,1);
TrainNum=round(row*(1-TestRatio));
TestNum=row-TrainNum;
TestIndices=zeros(TestNum,1);
for i=1:TestNum
    if TestRatio==0.2   
        step=1/TestRatio;
        TestIndices(i)=step*(i-1)+3;%5个中选第3个
    elseif TestRatio==0.25
        step=1/TestRatio;
        TestIndices(i)=step*(i-1)+2;%4个中选第2个     
    elseif TestRatio==1/3
        step=1/TestRatio;
        TestIndices(i)=step*(i-1)+2;%3个中选第2个 
    end 
end
spec_test = X(TestIndices,:);
target_test = Y(TestIndices,:);
A=ismember(X,spec_test);%相同元素行值全为1，不同元素行值全为0
B=A(:,1);
TrainIndices= ~B;
spec_train = X(TrainIndices,:);
target_train = Y(TrainIndices,:);
Data_Label.spec_train=spec_train;
Data_Label.target_train=target_train;
Data_Label.spec_test=spec_test;
Data_Label.target_test=target_test;
% step=3/TestRatio;
% con=fix(row/step);%取商
% rem=mod(row,step);%取余数
% for i=1:con+1
%     TestIndices((3*i-2))=step*(i-1)+3;%10个中选第3个
% %     if (3*i-2)<TestNum
% %        continue
% %     else 
% %         break
% %     end    
%     TestIndices((3*i-1))=step*(i-1)+6;%10个中选第6个
% %     if (3*i-2)<TestNum
% %        continue
% %     else 
% %         break
% %     end 
%     TestIndices(3*i)=step*(i-1)+9;%10个中选第6个
% %     if (3*i-2)<TestNum
% %        continue
% %     else 
% %         break
% %     end 
% end
end