function [Z96_train,Z96_test,A_mean,EIGENV,Z_train,Z_test] = PCA_decomp(train_features, test_features)
%3. Perform PCA analysis on the training set and project the test set onto the training set. 

%Shinjini Kundu (c) 2014
%Uses the PCA trick described in the paper Cootes TF, Taylor CJ, Cooper DH,
%Graham J, Active Shape Models - Their Training and Application. Computer
%Vision and Image Understanidng 1995 61(1): 38-59

%Inputs:       train_features   training matrix 
%              test_features    test matrix
%Outputs:      returns eigenvectors and feature matrix


NUM_TRAIN = size(train_features,1); 
NUM_TEST = size(test_features,1); 

A_mean = mean(train_features,1); %mean of each training feature across all training subjects
A = (train_features - repmat(A_mean,NUM_TRAIN,1)); %treating the data vectors as if they were in columns in a D * N fashion, although features matrix is N x D, where A represents the matrix X'
B = (test_features - repmat(A_mean,NUM_TEST,1)); 

T = (A*A')./NUM_TRAIN; 
[V,D] = eig(T); %compute orthonormal eigenvectors (V) and eigenvalues D such that Tv1 = dv1
%From eq. 41, A'*v1 is an eigenvector of S and has the same eigenvalue

V = fliplr(V); %order so that eigenvectors are in descreasing order
D = rot90(D,2); %order so that eigenvalues are in descending order

for i = 1:NUM_TRAIN
    eigenv = ((V(:,i)'*A)'./sqrt(D(i,i)*NUM_TRAIN)); %no need to make single and cause precision accuracy errors
    EIGENV(:,i) = eigenv; %column vectors of EIGENV
    %fprintf('Now computing eigenvector %d, which has norm %d and variance %d \n', i, norm(eigenv), D(i,i)); 
    Z_train(:,i) = A*eigenv; 
    Z_test(:,i) = B*eigenv; 
end

Z_train = real(Z_train);
Z_test = real(Z_test); 

%project only onto components that capture 96% of variance
variances = (cumsum(sort(diag(D),'descend'))./sum(diag(D))); 
for i = 1:NUM_TRAIN
    if floor(variances(i)*100)/100 >= 0.96
        cutoff = i; 
        break;
    end
end

%MATLAB does not give the same result between this code
%and the one below, to an approximation error. Will just use the code below for consistency with the
%other Run_PCA code
% Z96_train = A*EIGENV(:,1:cutoff);
% Z96_test = B*EIGENV(:,1:cutoff);

for i = 1:cutoff %change the number based on plot
    eigenv = EIGENV(:,i);
    Z96_train(:,i) = A*eigenv; 
    Z96_test(:,i) = B*eigenv; 
end



end

