function data = loaddata(num,rep)
%        读取光谱数据
%        通常情况下一个水果样本被采集多次，为此将被视为不同的样本进行建模
%        

SUM_absor = [];
SUM_spec = [];
data.meanabsor = [];
data.meanabsor = [];
for j=1:rep
    for i=1:num
        if i <10
        curfilename_A=['absor-00' int2str(i) '-' int2str(j) '.csv'];
        elseif i <100
            curfilename_A=['absor-0' int2str(i) '-' int2str(j) '.csv'];
        else
            curfilename_A=['absor-' int2str(i) '-' int2str(j) '.csv'];
        end
        curfile_A = csvread(curfilename_A);
        S_absor = curfile_A(:,2);
        SUM_absor =[SUM_absor,S_absor];
       
    end 
    
    for i=1:num
        if i <10
        curfilename_S=['spec-00' int2str(i) '-' int2str(j) '.csv'];
        elseif i <100
            curfilename_S=['spec-0' int2str(i) '-' int2str(j) '.csv'];
        else
            curfilename_S=['spec-' int2str(i) '-' int2str(j) '.csv'];
        end
        curfile = csvread(curfilename_S);  
        S_spec = curfile(:,2);
        SUM_spec =[SUM_spec,S_spec];
    end
    
end

data.wavelength = curfile(:,1)';
data.absor = SUM_absor';
data.spec = SUM_spec';


for k = 1:num
    data.meanabsor(k,:) = mean(data.absor(k:num:(num*rep),:),1);
    data.meanspec(k,:) = mean(data.spec(k:num:(num*rep),:),1);
end

