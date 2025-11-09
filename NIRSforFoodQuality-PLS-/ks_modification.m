function Data_Label=ks_modification(X,Y,Num) 
%  ks selects the samples XSelected which uniformly distributed in the exprimental data X's space  
%  Input   
%         X:the matrix of the sample spectra  
%         Y:the matrix of the sample measurement values
%         Num:the number of the sample spectra you want select
%  Output  
%         spec_train:the sample spectras was selected from the X, i.e. the X of the calibration set  
%         target_train:the sample measurement value was selected from the Y, i.e. the Y of the calibration set  
%         spec_test:the sample spectras remain int the X after select, i.e. the X of the prediction set 
%         target_test:the sample measurement value remain int the Y after select, i.e. the Y of the prediction set         
%  Raw Programmer: zhimin zhang @ central south university on oct 28 ,2007
%  Modification: Penghui Liu @ Zhejiang University on Oct 10, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% start of the kennard-stone step one  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    nR=size(X,1); % obtain the size of the X matrix  
    mDistance=zeros(nR,nR); %dim a matrix for the distance storage  
    vAllofSample=1:nR; 
    for i=1:nR-1  
        vRX=X(i,:); % 获取X的一行数据   
        for j=i+1:nR        
            vRX1=X(j,:); % 获得X中的另一行数据          
            mDistance(i,j)=norm(vRX-vRX1); % 计算欧氏距离                 
        end       
    end  
    [vMax,vIndexOfmDistance]=max(mDistance);  
    [~,nIndexofvMax]=max(vMax);   
    vSelectedSample(1)=nIndexofvMax;  
    vSelectedSample(2)=vIndexOfmDistance(nIndexofvMax);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% end of the kennard-stone step one  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% start of the kennard-stone step two  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    for i=3:Num  
        vNotSelectedSample=setdiff(vAllofSample,vSelectedSample);  
        vMinDistance=zeros(1,nR-i+1);       
        for j=1:(nR-i+1)  
            nIndexofNotSelected=vNotSelectedSample(j);  
            vDistanceNew = zeros(1,i-1); 
            for k=1:(i-1)  
                nIndexofSelected=vSelectedSample(k);  
                if(nIndexofSelected<=nIndexofNotSelected)  
                    vDistanceNew(k)=mDistance(nIndexofSelected,nIndexofNotSelected);  
                else  
                    vDistanceNew(k)=mDistance(nIndexofNotSelected,nIndexofSelected);      
                end                         
            end  
        vMinDistance(j)=min(vDistanceNew);  
        end     
        [~,nIndexofvMinDistance]=max(vMinDistance);  
        vSelectedSample(i)=vNotSelectedSample(nIndexofvMinDistance);  
    end  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %%%%% end of the kennard-stone step two  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %%%%% start of export the result  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    vSelectedRowIndex=vSelectedSample;  
    for i=1:length(vSelectedSample)    
        spec_train(i,:)=X(vSelectedRowIndex(i),:);  %训练集光谱数据
        target_train(i,:)=Y(vSelectedRowIndex(i),:);  %训练集性质数据
    end   
        vNotSelectedSample=setdiff(vAllofSample,vSelectedSample);  
    for i=1:length(vNotSelectedSample)      
        spec_test(i,:)=X(vNotSelectedSample(i),:);  %预测集光谱数据
        target_test(i,:)=Y(vNotSelectedSample(i),:);  %预测集性质数据
    end
Data_Label.spec_train=spec_train;
Data_Label.target_train=target_train;
Data_Label.spec_test=spec_test;
Data_Label.target_test=target_test;
end