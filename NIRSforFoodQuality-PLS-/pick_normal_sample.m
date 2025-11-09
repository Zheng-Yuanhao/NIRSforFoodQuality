tic
load('chucheng_EK336_600-950_468_20221210.mat')
II=[];
XX=[];
YY=[];
m=size(X,1);%样本数量
for i=1:m
    [pks,locs] = findpeaks(X(i,:));
    n=size(locs,2);%峰值数量
    for j=1:n
        if j<=n
            if locs(j)>210 & locs(j)<240    %EK336:210-230;EK343:210-230;
                x=X(i,:);
                y=Y(i,:);
                II=[II;i];
                XX=[XX;x];
                YY=[YY;y];
                i=[];
                x=[];
                y=[];
                j=n+1;
            else
                continue
            end
        end
        if j>n
            break
        end
    end
end
toc