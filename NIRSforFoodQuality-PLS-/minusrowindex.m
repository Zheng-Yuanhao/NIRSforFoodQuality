function [X_minus Y_minus]=minusrowindex(X,Y,minus_row_index)
    vSelectedRowIndex=minus_row_index';
    raw_row_index=linspace(1,size(X,1),size(X,1));
    vNotSelectedSample=setdiff(raw_row_index,minus_row_index);
    XX=zeros(size(vNotSelectedSample,2),size(X,2));
    YY=zeros(size(vNotSelectedSample,2),1);
    for i=1:size(vNotSelectedSample,2)     
            XX(i,:)=X(vNotSelectedSample(i),:);  %剩余光谱矩阵
            YY(i,:)=Y(vNotSelectedSample(i),:);  %剩余性质矩阵
    end
    X_minus=XX;
    Y_minus=YY;
end