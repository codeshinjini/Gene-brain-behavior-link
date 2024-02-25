function [T,V,D,Z,Z96,A_mean] = Run_PCA( features, matter,NUM_SUBJECTS,labels, keep )
%3. Perform PCA analysis and save the eigenvectors. The top eigenvector is
%called eigenvector_1 and corresponds to the largest variance in the matrix
%D. 

%Shinjini Kundu (c) 2014
%Uses the PCA trick described in the paper Cootes TF, Taylor CJ, Cooper DH,
%Graham J, Active Shape Models - Their Training and Application. Computer
%Vision and Image Understanidng 1995 61(1): 38-59

%Inputs:       features         data matrix (NxD)
%              NUM_SUBJECTS     N
%              matter           gm or wm
%              keep             1 for save, 0 for no save
%Outputs:      saves eigenvectors and feature matrix to folders

%clear X Y Z; 
%matter = 'medial';

A_mean = mean(features,1); %mean of each feature D across all subjects
A = (features - repmat(A_mean,NUM_SUBJECTS,1)); %treating the data vectors as if they were in columns in a D * N fashion, although features matrix is N x D, where A represents the matrix X'

% % try 
% %     [~,S,V] = svd(A,'econ');
% %     D = S.^2; 
% %     fprintf('SVD worked \n');
% %     Z = real(A*V); 
% %     
% %     if keep
% %         save(strcat(matter,'_feats'),'V','D','Z','A_mean');
% %     end
% %     
% % catch

    T = (A*A')./NUM_SUBJECTS; 
    [V,D] = eig(T); %compute orthonormal eigenvectors (V) and eigenvalues D such that Tv1 = dv1
    %from eq. 41, A'*v1 is an eigenvector of S and has the same eigenvalue

    V = fliplr(V); %order so that eigenvectors are in descreasing order
    D = rot90(D,2); %order so that eigenvalues are in descending order

    for i = 1:NUM_SUBJECTS
        eigenv = ((V(:,i)'*A)'./sqrt(D(i,i)*NUM_SUBJECTS)); %no need to make single and cause precision accuracy errors
        if keep
            save(strcat('eigenvector_',num2str(i)),'eigenv')
        end
        fprintf('Now computing eigenvector %d, which has norm %d and variance %d \n', i, norm(eigenv), D(i,i)); 
        Z(:,i) = A*eigenv; 
    end

    Z = real(Z); 

    if keep
        save(strcat(matter,'_feats'),'T','V','D','Z','A_mean');
    end
% end	

if keep
%plot the fraction of variance captured by the components
plot(cumsum(sort(diag(D),'descend'))./sum(diag(D)))
grid on; 
title('fraction of variance captured by component');
xlabel('principal component') 
ylabel('fraction of variance captured'); 
end

%project only onto components that capture 96% of variance
variances = (cumsum(sort(diag(D),'descend'))./sum(diag(D))); 
for i = 1:NUM_SUBJECTS
    if floor(variances(i)*100)/100 >= 0.96
        cutoff = i; 
        break;
    end
end

for i = 1:cutoff %change the number based on plot
   load(strcat('eigenvector_',num2str(i))); 
   Z96(:,i) = A*eigenv;
end

if keep
   save Z96 Z96 labels
end
%Z96 = [];

end

