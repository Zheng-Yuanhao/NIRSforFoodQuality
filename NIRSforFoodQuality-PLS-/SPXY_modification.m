function Data_Label=SPXY_modification(X,Y,Num) 
%  SPXY selects the samples XSelected which uniformly distributed in the exprimental data X's space and data Y's space  
%  Input   
%         X:the matrix of the sample spectra  
%         Y:the matrix of the sample measurement values
%         Num:the number of the sample spectra you want select,
%         eg:the number of sample set is 120,the test_size=0.3(7:3),and Num is 120*(1-0.3)=84. 
%  Output  
%         spec_train:the sample spectras was selected from the X, i.e. the X of the calibration set  
%         target_train:the sample measurement value was selected from the Y, i.e. the Y of the calibration set  
%         spec_test:the sample spectras remain int the X after select, i.e. the X of the prediction set 
%         target_test:the sample measurement value remain int the Y after select, i.e. the Y of the prediction set         
%  Programmer: zhimin zhang @ central south university on oct 28 ,2007
%  Modification: Penghui Liu @ Zhejiang University on Oct 10, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% start of the kennard-stone step one  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
dminmax = zeros(1,size(Y,1)); % Initializes the vector of minimum distances 
M = size(X,1); % Number of rows in X (samples) 
samples = 1:M; 
Dx = zeros(M,M); % Initializes the matrix of X distances 
Dy = zeros(M,M); % Initializes the matrix of Y distances 
for i=1:M-1 
    xa = X(i,:); 
    ya = Y(i); 
    for j = i+1:M 
      xb = X(j,:); 
      yb = Y(j); 
      Dx(i,j) = norm(xa - xb); 
      Dy(i,j) = norm(ya - yb); 
    end 
end 
Dxmax = max(max(Dx)); 
Dymax = max(max(Dy)); 
Dxy = Dx/Dxmax + Dy/Dymax; % Combines the distances in X and Y 
[maxD,index_row] = max(Dxy); % maxD = Row vector containing the largest element of each column in D 
[~,index_column] = max(maxD); % index_column = column corresponding to the largest element in matrix D 
m(1) = index_row(index_column); 
m(2) = index_column; 
dminmax(2) = Dxy(m(1),m(2)); 
for i = 3:Num
    % This routine determines the distances between each sample still available for selection and each of the samples already selected 
    pool = setdiff(samples,m); % pool = Samples still available for selection 
    dmin = zeros(1,M-i+1); % Initializes the vector of minimum distances between each sample in pool and the samples already selected 
    for j = 1:(M-i+1) % For each sample xa still available for selection 
        indexa = pool(j); % indexa = index of the j-th sample in pool (still available for selection) 
        d = zeros(1,i-1); % Initializes the vector of distances between the j-th sample in pool and the samples already selected 
        for k = 1:(i-1) % The distance with respect to each sample already selected is analyzed 
            indexb =  m(k); % indexb = index of the k-th sample already selected 
            if indexa < indexb 
                d(k) = Dxy(indexa,indexb); 
            else 
                d(k) = Dxy(indexb,indexa); 
            end 
        end 
        dmin(j) = min(d); 
    end 
    % The selected sample corresponds to the largest dmin 
    [dminmax(i),index] = max(dmin); 
    m(i) = pool(index); 
    m_complement = setdiff(linspace(1,M,M), m);
    spec_train = X(m, :); 
    target_train = Y(m); 
    spec_test = X(m_complement, :);
    target_test = Y(m_complement); 
end
Data_Label.spec_train=spec_train;
Data_Label.target_train=target_train;
Data_Label.spec_test=spec_test;
Data_Label.target_test=target_test;
end