%% results of the 10-fold cross validation 

addpath 'Analysis codes'
load demographic; load paper_indices; load Originals/original_volumes.mat;

%age and gender only 
Z = [(ages(inds==1)-mean(ages(inds==1)))./std(ages(inds==1)) (genders(inds==1)-mean(genders(inds==1)))./std(genders(inds==1)) [vol(inds==1)-mean(vol(inds==1))]'./std(vol(inds==1))]; 
L = group(inds==1)';


cross_validation_script(Z(:,1:2),L); %age and gender only

cross_validation_script(Z,L); %age, gender, total volumes

cross_validation_script(Z(:,3),L); %volumes only

% Visualize 
results = VisualizePLDA(Z96,labels,3,I0,A_mean); 
