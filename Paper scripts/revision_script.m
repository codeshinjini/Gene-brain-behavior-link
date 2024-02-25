%%% Revision materials script (c) 2024  %%
%%% This additional script helps generate some of the results requested by
%%% reviewers that are included in the supplementary materials.

%% Classification results when regressing out the effects of confounding
%variables of age and gender

matter = 'WM'; %can change depending on which result to run or GM

load demographic; %load demographic variables
load paper_indices; %indices to select which subjects after exclusion criteria

Z = [ages(inds==1) genders(inds==1)]; %confounding variables
if isequal(matter,'GM')
    clear features; 
    load TBM_gray_down/PCA_TBM/features.mat; %load features pathname
elseif isequal(matter,'WM')
    clear features; 
    load TBM_white_down/PCA_TBM/features.mat; %load features pathname
end

%according to Section S4 in the Supplementary materials

V = features - Z*((Z'*Z)\Z')*features; %residual of data matrix removing confounds

cd Supplementary; 
save('features_no_age_gender', 'label', 'V',"-v7.3");


%% run PCA
addpath 'Analysis codes'; 

[T,V,D,Z,Z96,A_mean] = Run_PCA( V, matter,206,label, 1 ); 

%% run PLDA

addpath 'Analysis codes'; 
[accuracy, kappa, sensitivity, specificity, PLDAprojected, w] = PLDA_knn(Z96, label+2, 1, 3, 'brain');

%% anatomic total image, normalized and smoothed with sigma = 1

load demographic; load paper_indices; sigma = 1;  

[Xt,Yt,Zt]=meshgrid(-3*sigma:3*sigma,-3*sigma:3*sigma,-3*sigma:3*sigma);
phi = gaussian_bf(Xt,Yt,Zt,sigma); %normalized gaussian kernel in 3D
%

indexes = find(inds==1);  

for i=1:numel(indexes) 
    load(strcat('registered/affine_down_skullstripT1_',num2str(indexes(i)),'.mat')); %segmented, affine registered whole brain images
    D = mat2gray(convn(skull_stripped,phi,'same')); %normalization
    features(i,:) = D(:)'; 
    i
end

%% Classification results with only affine registered normalized, smoothed images

sigma = 1;  

[Xt,Yt,Zt]=meshgrid(-3*sigma:3*sigma,-3*sigma:3*sigma,-3*sigma:3*sigma);
phi = gaussian_bf(Xt,Yt,Zt,sigma); %normalized gaussian kernel in 3D
%

% Anatomic image (smoothed and zeropadded, normalized in same manner as original) concatenate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matter = 'WM'; %or GM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indexes = find(inds==1); 
cd(matter); %or WM or GM
for i = 1:numel(indexes)
    if isequal(matter,'WM')
        load(strcat('c2affine_downT1_',num2str(indexes(i)),'.mat')); %resampled, affine registered, white matter maps
        zp = [15 25 15];
    elseif isequal(matter,'GM')
        load(strcat('c1affine_downT1_',num2str(indexes(i)),'.mat')); %resampled, affine registered, gray matter maps
        zp = [25 35 25];
    end
    D = mat2gray(convn(padarray(Data,zp),phi,'same'));
    features(i,:) = [D(:)']; 
    i
end

%save features features label -v7.3; 
%[T,V,D,Z,Z96,A_mean] = Run_PCA( features, 'Image',size(features,1),label, 1 );

%% Separate gray matter and white matter original volumes when doing age, gender, volume classification 

cd Paper scripts;
addpath 'Analysis codes'
load demographic; load paper_indices; load Originals/original_volumes.mat;

VOL = gm_vol; %which volume, choose vol gm_vol wm_vol

%age, gender, volume 
Z = [(ages(inds==1)-mean(ages(inds==1)))./std(ages(inds==1)) (genders(inds==1)-mean(genders(inds==1)))./std(genders(inds==1)) [VOL(inds==1)-mean(VOL(inds==1))]'./std(VOL(inds==1))]; 
L = group(inds==1)';

cross_validation_script(Z,L); %age, gender, total volumes

%% Calculate the classification tables

cross_validation_TBM('WM_orig'); 

%%%%%% To visualize image-domain maps %%%%%%%%%%%%


%% Average direction over 1000 iterations

direction1 = 0; direction2 = 0; 
A_mean = 0; 
T = 1000; 

for i = 1:T
    load(strcat('iter',num2str(i))); 
    direction1 = direction1 + w.dir1; 
    direction2 = direction2 + w.dir2; 
    A_mean = A_mean + w.MEAN; 
    i
    clear w;
end
direction1 = direction1/T; 
direction2 = direction2/T; 
A_mean = A_mean/T; 

LABELS = []; PROJECTED = []; 
for i = 1:T
        load(strcat('iter',num2str(i))); 
    for j = 1:10
        LABELS = [LABELS; w.train_labels{j}']; 
        PROJECTED = [PROJECTED; w.projected{j}];
    end
    i
end

save averaged_w direction1 direction2 A_mean LABELS PROJECTED; 

%% Slices to compare to image S6 in supplementary materials, only visualizing across direction 1
k = 3; 
slices = [55 81 90];

sigma1 = std(PROJECTED(:,1)); 
n1 = [-3 -1.5 0 1.5 3];

for j = 1:k
    pi(j) = sum(LABELS==j)/numel(LABELS);  
    mu_hat(j,:) = sum(PROJECTED((LABELS==j),:))./sum(LABELS==j);
end

for i = 1:numel(n1)
    fprintf('Now producing another brain image n1 = %d \n', n1(i)); 
    direction  = A_mean' + n1(i)*sigma1*direction1; %only visualizing across direction 1

    IMAGE{i} = imrotate(reshape(direction,[155 196 154]),90); 

    for m = 1:k
       distances(m) = norm([n1(i)*sigma1 0]-mu_hat(m,:),2)^2-log(pi(m));
    end
    [~,p] = min(distances); 
    membership(i) = p; 
end

[B,T] = Gen_Stack(IMAGE); 
results.T = T; %images are in the same orer as 
viewdicom(T,0); imcontrast; %the placement of the membership matrix is in the same placement as the images
%y axis corresponds to direction 2


BIG = [T(:,:,slices(3)); T(:,:,slices(2)); T(:,:,slices(1))]; 
imagesc(BIG); colorbar; colormap(jet)
ylabel('Axial slices'); 
ytick = size(BIG,1)/3; 
set(gca,'YTick',[ytick/2:ytick:3*ytick]); 
set(gca,'YTickLabel',{num2str(slices(3)),num2str(slices(2)),num2str(slices(1))});
set(gca,'XTick',[size(BIG,2)/10:size(BIG,2)/5:size(BIG,2)])
set(gca,'XTickLabel',[{'-3\sigma_1', '-1.5\sigma_1', '0', '1.5\sigma_1', '3\sigma_1'}]);
xlabel('Discriminant direction 1'); 
set(gca,'FontName','Times New Roman','FontSize',18);

imcontrast; 

