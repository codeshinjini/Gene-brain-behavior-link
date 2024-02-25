%% Shinjini Kundu (c) 2022
%% Script that helps generate the results in the paper

%% 1. Demographics
load demographic; %load the demographic variables  
load paper_indices; %indices in the filename order that are included in the paper
load original_volumes; %volumes of the original images


GROUP = 1; %choose -1, 0, or 1 for the group membership

indexes = find(inds==1); 

cohort_label = group(indexes); 
cohort_ages = ages(indexes); 

cohort_genders = genders(indexes); 
cohort_FSIQ = FSIQ(indexes); 
cohort_VIQ = VIQ(indexes); 
cohort_NVIQ = NVIQ(indexes); 
cohort_SRS = SRS(indexes); 
cohort_volume = vol(indexes); 

fprintf('Group total is %d \n', sum(cohort_label==GROUP));
fprintf('Mean age is %0.1f with std dev %0.1f and p value %0.2f \n', mean(cohort_ages(cohort_label==GROUP),'omitnan'), std(cohort_ages(cohort_label==GROUP),'omitnan'),anova1(cohort_ages,cohort_label,'off')); 
fprintf('Males %d and females %d and p value %0.2f \n', sum(cohort_genders(cohort_label==GROUP)==1), sum(cohort_genders(cohort_label==GROUP)==2),anova1(cohort_genders,cohort_label,'off')); 
fprintf('Mean volume is %.2f with std dev %0.2f and p value %.2f \n', mean(cohort_volume((cohort_label==GROUP))/10^6,'omitnan'), std(cohort_volume((cohort_label==GROUP))/10^6,'omitnan'),anova1(cohort_volume,cohort_label,'off'));

fprintf('Mean FSIQ is %.2f with std dev %.2f and p value %.2f \n', mean(cohort_FSIQ(cohort_label==GROUP),'omitnan'), std(cohort_FSIQ(cohort_label==GROUP),'omitnan'),anova1(cohort_FSIQ,cohort_label,'off')); 
fprintf('Mean VIQ is %.2f with std dev %0.2f and p value %0.2f \n', mean(cohort_VIQ(cohort_label==GROUP),'omitnan'), std(cohort_VIQ(cohort_label==GROUP),'omitnan'),anova1(cohort_VIQ,cohort_label,'off')); 
fprintf('Mean NVIQ is %0.2f with std dev %0.2f and p value %0.2f \n', mean(cohort_NVIQ(cohort_label==GROUP),'omitnan'), std(cohort_NVIQ(cohort_label==GROUP),'omitnan'), anova1(cohort_NVIQ,cohort_label,'off')); 

fprintf('Mean SRS is %0.2f with std dev %0.2f and p value %0.2f \n', mean(cohort_SRS(cohort_label==GROUP),'omitnan'), std(cohort_SRS(cohort_label==GROUP),'omitnan'), anova1(cohort_SRS,cohort_label,'off')); 


find(isnan(cohort_ages(cohort_label==GROUP)))
find(isnan(cohort_volume(cohort_label==GROUP)))
find(isnan(cohort_FSIQ(cohort_label==GROUP)))
find(isnan(cohort_VIQ(cohort_label==GROUP)))
find(isnan(cohort_NVIQ(cohort_label==GROUP)))
find(isnan(cohort_SRS(cohort_label==GROUP)))


%% 2. Calculate numbers with diagnoses

%duplications and deletions occupy first 1-88 spots on cohort_label
DUPDEL = 88; 

cADHD = ADHD(indexes); cANXIETY = ANXIETY(indexes); cARTICULATION = ARTICULATION(indexes); 
cBEHAVIOR = BEHAVIOR(indexes); cASD = ASD(indexes); cCOORDINATION = COORDINATION(indexes); 
cMOOD = MOOD(indexes); cMR = MR(indexes); cSTEREOTYPE = STEREOTYPE(indexes);
cTIC = TIC(indexes); cENURESIS = ENURESIS(indexes); cLANGUAGE = LANGUAGE(indexes); cLEARNING = LEARNING(indexes);
fprintf('Total %d with ADHD, or %0.2f percent \n', sum(cADHD(1:DUPDEL)==1), sum(cADHD(1:DUPDEL)==1)/DUPDEL*100);
fprintf('Total %d with ANXIETY, or %0.2f percent \n', sum(cANXIETY(1:DUPDEL)==1), sum(cANXIETY(1:DUPDEL)==1)/DUPDEL*100); 
fprintf('Total %d with ARTICULATION, or %0.2f percent \n', sum(cARTICULATION(1:DUPDEL)==1), sum(cARTICULATION(1:DUPDEL)==1)/DUPDEL*100); 
fprintf('Total %d with BEHAVIOR, or %0.2f percent \n', sum(cBEHAVIOR(1:DUPDEL)==1), sum(cBEHAVIOR(1:DUPDEL)==1)/DUPDEL*100); 
fprintf('Total %d with ASD, or %0.2f percent \n', sum(cASD(1:DUPDEL)==1), sum(cASD(1:DUPDEL)==1)/DUPDEL*100); 
fprintf('Total %d with COORDINATION, or %0.2f percent \n', sum(cCOORDINATION(1:DUPDEL)==1), sum(cCOORDINATION(1:DUPDEL)==1)/DUPDEL*100); 
fprintf('Total %d with TIC, or %0.2f percent \n', sum(cTIC(1:DUPDEL)==1), sum(cTIC(1:DUPDEL)==1)/DUPDEL*100); 
fprintf('Total %d with ENURESIS, or %0.2f percent \n', sum(cENURESIS(1:DUPDEL)==1), sum(cENURESIS(1:DUPDEL)==1)/DUPDEL*100); 
fprintf('Total %d with LANGUAGE, or %0.2f percent \n', sum(cLANGUAGE(1:DUPDEL)==1), sum(cLANGUAGE(1:DUPDEL)==1)/DUPDEL*100); 
fprintf('Total %d with LEARNING, or %0.2f percent \n', sum(cLEARNING(1:DUPDEL)==1), sum(cLEARNING(1:DUPDEL)==1)/DUPDEL*100); 
fprintf('Total %d with MOOD, or %0.2f percent \n', sum(cMOOD(1:DUPDEL)==1), sum(cMOOD(1:DUPDEL)==1)/DUPDEL*100); 
fprintf('Total %d with MR, or %0.2f percent \n', sum(cMR(1:DUPDEL)==1), sum(cMR(1:DUPDEL)==1)/DUPDEL*100); 
fprintf('Total %d with STEREOTYPE, or %0.2f percent \n', sum(cSTEREOTYPE(1:DUPDEL)==1), sum(cSTEREOTYPE(1:DUPDEL)==1)/DUPDEL*100); 

meanDIAGS = cADHD + cANXIETY + cARTICULATION + cBEHAVIOR + cASD + cCOORDINATION + cMOOD + cMR + cSTEREOTYPE +cTIC + cENURESIS + cLANGUAGE + cLEARNING;


fprintf('Total missing %d \n', sum(isnan(meanDIAGS(1:DUPDEL))));
[mean(meanDIAGS(~isnan(meanDIAGS(1:DUPDEL)))) std(meanDIAGS(~isnan(meanDIAGS(1:DUPDEL))))]


%% 4. Run TBM on independent template

%parameters to run 1 brain successfully from initialization of identity
%numScales = 4; lambda = 3; stepsize = 10^-7, 10^-5, 10^-1, cutoff is
%4*10^-4 and then 5*10^-3

addpath TBM_codes; 

load paper_indices.mat

%white matter processing
cd WM; load I0_white_down_new.mat; paddedI0 = padarray(I0_white_down_new,[15 25 15]); %for gray matter, padding is [25 35 25]

%execute for each brain image:  
new_results = multVOT_init(paddedI0,paddedI1,results.f1,results.f2,results.f3);
        
  

paddedI1 = padarray(Data,[25 35 25]); %for gray matter  


%% 5. Concatenate all features

load paper_indices; 
load demographic; 
addpath 'Analysis codes'
load T1_1_results %results of TBM processed subject 1
I0_template = results.I0; 
[X,Y,Z] = meshgrid(1:size(results.I0,2),1:size(results.I0,1),1:size(results.I0,3)); 
indexes = find(inds==1);
counter = 1;
for i = 1:numel(indexes)
    load(strcat('T1_',num2str(indexes(i)),'_results'));
    features(counter,:) = single([results.f1(:)-X(:); results.f2(:)-Y(:); results.f3(:)-Z(:)]'); 
    clear results;
    label(counter) = group(indexes(i));
    counter = counter + 1;
    i
end

save features features label -v7.3

%% 6. Run PCA

[T,V,D,Z,Z96,A_mean] = Run_PCA( features, 'WM',206,label, 1 ); 

%% 
%generate PCA plots paper
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matter = 'WM'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%change directory to where the TBM processed files are
load PCA_Image_raw/Image_feats; 
plot(cumsum(sort(diag(D),'descend'))./sum(diag(D)),'LineWidth',2)
hold on; clear D; 
if isequal(matter,'WM')
    load PCA_TBM/WM_feats; 
elseif isequal(matter,'GM')
    load PCA_TBM/GM_feats; 
end
plot(cumsum(sort(diag(D),'descend'))./sum(diag(D)),'LineWidth',2)
legend('image domain','transport domain')
grid on; 
if isequal(matter,'WM')
    title('White matter')
elseif isequal(matter,'GM')
    title('Gray matter')
    ylim([0 1.1]);
end

xlim([0 206])
set(gca,'FontSize',18); 
set(gca,'fontname','Times New Roman');
ylabel('Fraction of variance'); 
xlabel('Number of components')

%% 7 . PLDA plots

addpath 'Analysis codes'; 
[accuracy, kappa, sensitivity, specificity, PLDAprojected, w] = PLDA_knn(Z96, labels+2, 1, 3, 'brain');

%can use the PLDAprojected scores for behavior correlation

%% For cross-validation 

classification_script_paper; 

%% To create histograms of the original brain volumes
load Originals/original_volumes.mat;
load demographic; load paper_indices; L = group(inds==1); V = vol(inds==1)'; 

for i = -1:1
    i
    h = histfit(V(L==i)); 
    CurveX{i+2} = h(2).XData;
    CurveY{i+2} = h(2).YData;
end
color(1,:) = [0.5 0.5 0.5]; color(2,:) = [0 0.4 0.4]; color(3,:) = [1 0.5 0.1]; 

figure; hold on
for i = 1:3
    plot(CurveX{i}/10^6,CurveY{i},'LineWidth',2,'color',color(i,:)); 
end

set(gca,'FontSize',18); 
set(gca,'fontname','Times New Roman');grid on;
ylabel('Number of individuals'); xlabel('Brain tissue volume (L)'); 
legend('Deletion Carriers','Controls','Duplication Carriers'); 


%% 8. 10-fold CV repeated 1000x
 
%generate folds
load paper_indices; 
indexes = find(inds==1); 

nFolds = 10;
c = cvpartition(labels,'KFold',nFolds,'Stratify',true); %initial partition
C{1} = c;
rng(1); %set seed for random number generator
for i = 1:999
    C{i+1} = repartition(c); 
    i
end

%confirming that the partitions are disjoint
for i = 1:1000
    for j = 1:1000
        if i==j
            e2(i,j) = 0; 
        else
            e2(i,j) = isequal(C{i},C{j}); 
        end
    end
end

save CVpartitions C; %saved under classification folder 


%% PLDA and CV
%go to cross_validation_script in the main folder 
  
%cross-validation 
cross_validation_script('WM');
cross_validation_TBM('GM');

%% %% Harvard-Oxford Atlas on displacement maps to determine areas of greatest shift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%matter = 'WM'; 
n = 3; %0 3 or -3 %standard deviations from the mean
WM_REG = [1 3 8 12 14]; %subcortical atlas values for white matter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:2
    if k==1
        matter = 'GM'; 
    elseif k==2
        matter = 'WM'; 
    end
%load the displacement fields
cd Visualization; %where visualizations are saved
load(strcat(matter,'_TBM_direction.mat')); 

%calculate the scores for standard deviations
sigma1 = std(PROJECTED(:,1)); 
sigma2 = std(PROJECTED(:,2)); 

%load the corresponding template image, in the same space as the template
if isequal(matter,'GM')
    load TBM_gray_down/T1_1_results.mat; %load TBM processed gray matter subject 1
    I0_trunc = imrotate(results.I0(26:end-25,36:end-35,26:end-25),90);
elseif isequal(matter,'WM')
    load TBM_white_down/T1_1_results.mat; %load TBM processed white matter subject 1
    I0_trunc = imrotate(results.I0(16:end-15,26:end-25,16:end-15),90);
end
[M,N,K] =size(results.I0);
[X,Y,Z] = meshgrid(1:N,1:M,1:K); 



%calculate displacement maps for both directions 1 and 2
null_mean = 0; %0 to exclude it, 1 to include it
for i = 1:2
     if i==1
         direction = A_mean'*null_mean +n*sigma1*direction1;
     elseif i==2
         direction = A_mean'*null_mean +n*sigma2*direction2;
     end
     sz = size(direction,1)/3;
     u = reshape(direction(1:sz),M,N,K);
     v = reshape(direction(sz+1:2*sz),M,N,K); 
     z = reshape(direction(2*sz+1:3*sz),M,N,K); 
    
    f = double(X + u); 
    g = double(Y + v); 
    h = double(Z + z); 

    [dfdx,dfdy,dfdz] = gradient(f); %compute Jacobian map
    [dgdx,dgdy,dgdz] = gradient(g); 
    [dhdx,dhdy,dhdz] = gradient(h);

    %And Jacobian determinant |Du|
    D = (dfdx.*dgdy.*dhdz + dfdy.*dgdz.*dhdx + dfdz.*dgdx.*dhdy - dfdx.*dgdz.*dhdy - dfdy.*dgdx.*dhdz - dgdy.*dhdx.*dfdz); %determinant

     if i==1
         if isequal(matter,'GM')
            jacobian1 = imrotate(D(26:end-25,36:end-35,26:end-25),90); 
         elseif isequal(matter,'WM')
            jacobian1 = imrotate(D(16:end-15,26:end-25,16:end-15),90); 
         end
         clear D;
     elseif i==2
         if isequal(matter,'GM')
            jacobian2 = imrotate(D(26:end-25,36:end-35,26:end-25),90); 
         elseif isequal(matter,'WM')
            jacobian2 = imrotate(D(16:end-15,26:end-25,16:end-15),90); 
         end
         clear D;
     end
 end
 
%load atlases 
load atlases; %load the cortical and subcortical atlases
CORTICAL = round(cortical_atlas); SUBCORTICAL = round(subcortical_atlas); %round atlases to whole numbers
mask_cortical = (cortical_atlas~=0); %mask for gray matter images only 
mask_subcortical =(subcortical_atlas~=0); %mask for gray matter and white matter

%calculate normalized displacement maps
if isequal(matter,'WM')
    dmap1 = mask_subcortical.*jacobian1; 
    dmap2 = mask_subcortical.*jacobian2; 
 
    Zdmap1 = (dmap1-mean(dmap1(dmap1~=0)))./std(dmap1(dmap1~=0)); %zscore map
    Zdmap2 = (dmap2-mean(dmap2(dmap2~=0)))./std(dmap2(dmap2~=0)); %zscore map

    %subcortical atlas
    for i = 1:5
        mask = medfilt3(SUBCORTICAL==WM_REG(i)); %median filter to null edge signal
        SCD1(WM_REG(i)) = sum(Zdmap1(:).*mask(:))./nnz(mask); %vector of normalized displacements for direction 1
        SCD2(WM_REG(i)) = sum(Zdmap2(:).*mask(:))./nnz(mask); %vector of normalized displacmeents for direction 2
        clear mask;
    end
elseif isequal(matter,'GM')
    dmap1 = mask_cortical.*jacobian1; 
    dmap2 = mask_cortical.*jacobian2; 
     
    Zdmap1 = (dmap1-mean(dmap1(dmap1~=0)))./std(dmap1(dmap1~=0)); %zscore map
    Zdmap2 = (dmap2-mean(dmap2(dmap2~=0)))./std(dmap2(dmap2~=0)); %zscore map

    %cortical atlas
    for i = 1:48
        mask = medfilt3(CORTICAL==i); %median filter to null edge signal
        CD1(i) = sum(sum(sum(Zdmap1.*mask)))./nnz(mask); %vector of normalized displacements for direction 1
        CD2(i) = sum(sum(sum(Zdmap2.*mask)))./nnz(mask); %vector of normalized displacmeents for direction 2
        clear mask; 
    end

    dmap1 = mask_subcortical.*jacobian1; 
    dmap2 = mask_subcortical.*jacobian2; 

    Zdmap1 = (dmap1-mean(dmap1(dmap1~=0)))./std(dmap1(dmap1~=0)); %zscore map
    Zdmap2 = (dmap2-mean(dmap2(dmap2~=0)))./std(dmap2(dmap2~=0)); %zscore map
     %subcortical atlas
    for i = 1:21
        if ~any(WM_REG==i)
            mask = medfilt3(SUBCORTICAL==i); %median filter to null edge signal
            SCD1(i) = sum(sum(sum(Zdmap1.*mask)))./nnz(mask); %vector of normalized displacements for direction 1
            SCD2(i) = sum(sum(sum(Zdmap2.*mask)))./nnz(mask); %vector of normalized displacmeents for direction 2
            clear mask; 
        end
    end

end

end

round(CD2'*100)/100;
[M,D] = sort(CD2,'descend'); 
[M(1:3)' D(1:3)']
[M(end-2:end)' D(end-2:end)']

% round(SCD2'*100)/100;
% [M,D] = sort(SCD2,'descend'); 
% [M' D']

%% Regression analysis: Appendix volume vs. tissue distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matter = 'GM';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath 'Analysis codes';

indexes = find(inds==1); exclude = find(isnan(FSIQ(indexes))); indexes(exclude) = [];

if isequal(matter,'GM')
    load TBM_gray_down/PCA_TBM/Z96;
    load TBM_gray_down/PCA_TBM/GM_feats;
    load TBM_gray_down/T1_1_results;
    VOLUMES = gm_vol(indexes);
elseif isequal(matter,'WM')
    load TBM_white_down/PCA_TBM/Z96; 
    load TBM_white_down/PCA_TBM/WM_feats;
    load TBM_white_down/T1_1_results;
    VOLUMES = wm_vol(indexes);
end

features = Z96; features(exclude,:) = [];

results = Run_Regression(features,VOLUMES'/10^6,1,results.I0,A_mean,'TBM',0,[ages(indexes) genders(indexes) group(indexes) FSIQ(indexes)]); 
xlabel('Brain tissue volume (L)'); 
if isequal(matter,'WM')
    title('white matter volume vs. distribution')
elseif isequal(matter,'GM')
    title('gray matter volume vs. distribution')
end

%text box label (optional)
dim = [.2 .5 .3 .3];
if results.p==0
    str = {strcat('R = ',num2str(round(results.CC_train*10^2)/10^2)), strcat('p <0.001 ')}; %= ',num2str(results.p))};
else
    str = {strcat('R = ',num2str(round(results.CC_train*10^2)/10^2)), strcat('p = ',num2str(round(results.p*10^2)/10^2))}; %= ',num2str(results.p))};
end
annotation('textbox',dim,'String',str,'position',[0.73 0.25,0.1,0.1],'FitBoxToText','on');
%white matter and FSIQ, VIQ, NVIQ is positive, correcting for age and genders

%2. Correlation between gray and white matter projection

%% Determinant of Jacobian
%% this is the code used to generate the jacobian maps in the jacobians.mat file
%direction 1 only 
B = []; %i = 2; 
for i = 1:2
    [M,N,K] = size(results.I0); 

%calculate the scores for standard deviations
sigma1 = std(PROJECTED(:,1)); 
sigma2 = std(PROJECTED(:,2)); 
if i==1
n1 = -3; %interested in examining the shifts 1.5 standard deviations from the mean
elseif i==2
    n1=3; 
end
n2 = 0; 
direction = A_mean' +n1*sigma1*direction1 + n2*sigma2*direction2;


 sz = size(direction,1)/3;
 u = reshape(direction(1:sz),M,N,K);
 v = reshape(direction(sz+1:2*sz),M,N,K); 
 z = reshape(direction(2*sz+1:3*sz),M,N,K); 
 [X,Y,Z] = meshgrid(1:N,1:M,1:K); 

 f = double(X + u); 
 g = double(Y + v); 
 h = double(Z + z); 

 [dfdx,dfdy,dfdz] = gradient(f); %compute Jacobian map
 [dgdx,dgdy,dgdz] = gradient(g); 
 [dhdx,dhdy,dhdz] = gradient(h);

 %And Jacobian determinant |Du|
 D = (dfdx.*dgdy.*dhdz + dfdy.*dgdz.*dhdx + dfdz.*dgdx.*dhdy - dfdx.*dgdz.*dhdy - dfdy.*dgdx.*dhdz - dgdy.*dhdx.*dfdz); %determinant

% %%gray matter mask
CONST = 1.5;
 mask = results.I0>CONST*min(results.I0); 

 masked = mask.*D; 
 Zmasked = (masked-mean(masked(masked~=0)))./std(masked(masked~=0)); %zscore map
B = [B imrotate(Zmasked,90)];

if i==1
    del_D = D; 
elseif i==2
    dup_D = D; 
end

end
 close all; 
 viewdicom(B,0)

 %calculate correlation between 
 %since pvalue cannot be "zero" will report "0" as p<0.001
 [rho, pval] = corr(dup_D(:),del_D(:))

 %% this is to rotate the maps after the jacobian above is calculated
mode = 'WM';
map = del_D; 

%%%%%%%%%%%%%code is same after that%%%%%%
for i = 1:104
    if isequal(mode,'WM')
        map_rotated = imrotate(map(16:end-15,26:end-25,16:end-15),90); 
    end
    if isequal(mode,'GM')
        map_rotated = imrotate(map(26:end-25,36:end-35,26:end-25),90); 
    end
end



%% Visualize whole image panel
n1 = [-3 -1.5 0 1.5 3];n2 = [-3 -1.5 0 1.5 3];
[B,T] = Gen_Stack(image(1:5,1:5)); 
%viewdicom(T,0); 

slice = 55;
imagesc(T(:,:,slice)); %colorbar; 
set(gca,'FontSize',18); 
set(gca,'FontName','Times New Roman'); 
set(gcf, 'Position', [0 0 750 750]); 
xlabel('Discriminant direction 1'); 
ylabel('Discriminant direction 2');
colormap(jet); colorbar;
xticks(size(T,2)/[numel(n1)*2]:size(T,2)/numel(n1):size(T,2))
yticks(size(T,1)/[numel(n2)*2]:size(T,1)/numel(n2):size(T,1))
arr1 = []; arr2 = [];
for i = 1:numel(n1)
    if n1(i)==0
            arr1 = [arr1, '0' ];
    else
    arr1 = [arr1, {strcat(num2str(n1(i)),'\sigma_1')} ];
    end
end
for i = 1:numel(n2)
    if n2(i)==0
        arr2 = [arr2, '0'];
    else
    arr2 = [arr2, {strcat(num2str(-n2(i)),'\sigma_2')}]; 
    end
end
xticklabels(arr1); yticklabels(arr2);
colormap(jet)
%xticklabels([{strcat(num2str(n1(1)),'\sigma_1')},{strcat(num2str(n1(2)),'\sigma_1')},{strcat(num2str(n1(3)))}, {strcat(num2str(n1(4)),'\sigma_1')}, {strcat(num2str(n1(5)),'\sigma_1')}]); 
%yticklabels([{strcat(num2str(n2(5)),'\sigma_2')},{strcat(num2str(n2(4)),'\sigma_2')},{strcat(num2str(n2(3)))}, {strcat(num2str(n2(2)),'\sigma_2')}, {strcat(num2str(n2(1)),'\sigma_2')}]); 

%% Visualize individual axial slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matter = 'WM'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(matter,'WM')
    load Visualization/WM_TBM_visualization.mat; %where visualization files are saved
elseif isequal(matter,'GM')
    load Visualization/GM_TBM_visualization.mat; %where visualization files are saved
end

[B,T] = Gen_Stack(image(3,:)); %grabbing images only from direction 1 
%viewdicom(T,0); 

if isequal(matter,'WM')
    slices = [85,73,63,53,47,39];
elseif isequal(matter,'GM')
    slices = [90,81,55];
end

BIG = []; 
yl = [];
for i = 1:numel(slices)
    BIG = [BIG; T(:,:,slices(i))]; 
    yl = [yl {strcat(num2str(slices(i)))}]; 
end
imagesc(BIG); colormap(jet); colorbar; 
set(gca,'FontSize',18); 
set(gca,'FontName','Times New Roman'); 
xlabel('Discriminant direction 1'); 
ylabel('Axial slice');
colormap(jet)
xticks(size(T,2)/10:size(T,2)/5:size(T,2))
yticks(size(T,1)/2:size(T,1):(numel(slices))*size(T,1))
xticklabels([{strcat(num2str(n1(1)),'\sigma_1')},{strcat(num2str(n1(2)),'\sigma_1')},{strcat(num2str(n1(3)))}, {strcat(num2str(n1(4)),'\sigma_1')}, {strcat(num2str(n1(5)),'\sigma_1')}]); 
yticklabels(yl); 
set(gcf, 'Position', [0 0 750 750]); 


%% A. Multiple linear regression: behavioral correlations (training set only, with p-value)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matter = 'GM'; %can change
VARIABLE = NVIQ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(matter,'GM')
    load TBM_gray_down/PCA_TBM/Z96; 
end
if isequal(matter,'WM')
    load TBM_white_down/PCA_TBM/Z96; 
end
addpath 'Analysis codes'; 
[accuracy, kappa, sensitivity, specificity, PLDAprojected, w] = PLDA_knn(Z96, labels+2, 1, 3, 'brain');
close all;

% load the demographics
load('demographic'); 
load('paper_indices'); 
 
var = VARIABLE(inds==1); %articulation, coordination, enuresis, language, learning
projected_scores = PLDAprojected; 
%age, gender and volumes already proven to not be enough to differentiate
%the groups, not including in the model

%eliminate nans
toremove = find(isnan(var)); %v(toremove) = [];
var(toremove) = []; projected_scores(toremove,:) = []; 
%ag(toremove) = []; ge(toremove) = [];

%now compute correlations
[b,bint,r,rint,stats] = regress(var,[ones(size(projected_scores,1),1) projected_scores]); %syntax according to documented 
R_squared = stats(1); 
coefficients = bint; 

%calculate p-value
rng(1); count = 0; 
for i = 1:1000
    r = randsample(numel(var),numel(var)); 
    [b,bint,r,rint,stats] = regress(var(r),[ones(size(projected_scores,1),1) projected_scores]); %syntax according to documented 
    if stats(1) >= R_squared
        count = count  +1;
    end
end
p = count./1000;

fprintf('Score was missing for %1.1f subjects \n',numel(toremove));
fprintf('R^2 is %1.2f with p value %1.2f \n',R_squared,p)
coefficients


%% mulitple linear regression: combining white matter and gray matter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VARIABLE = NVIQ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:2
    if i==1
    load TBM_gray_down/PCA_TBM/Z96; 
    elseif i==2
    load TBM_white_down/PCA_TBM/Z96; 
    end
addpath 'Analysis codes'; 
[accuracy, kappa, sensitivity, specificity, PLDAprojected, w] = PLDA_knn(Z96, labels+2, 1, 3, 'brain');
close all;
if i==1
    GMproj = PLDAprojected; 
elseif i==2
    WMproj= PLDAprojected; 
end
end
% load the demographics
load('demographic'); 
load('paper_indices'); 
 
var = VARIABLE(inds==1); %articulation, coordination, enuresis, language, learning
projected_scores_GM = GMproj; projected_scores_WM = WMproj; 
%age, gender and volumes already proven to not be enough to differentiate
%the groups, not including in the model

%eliminate nans
toremove = find(isnan(var)); %v(toremove) = [];
var(toremove) = []; projected_scores_GM(toremove,:) = []; projected_scores_WM(toremove,:) = []; 
%ag(toremove) = []; ge(toremove) = [];

%now compute correlations
[b,bint,r,rint,stats] = regress(var,[ones(size(projected_scores,1),1) projected_scores_GM projected_scores_WM]); %syntax according to documented 
R_squared = stats(1); 
coefficients = bint; 

%calculate p-value
rng(1); count = 0; 
for i = 1:1000
    r = randsample(numel(var),numel(var)); 
    [b,bint,r,rint,stats] = regress(var(r),[ones(size(projected_scores,1),1) projected_scores_GM projected_scores_WM]); %syntax according to documented 
    if stats(1) >= R_squared
        count = count  +1;
    end
end
p = count./1000;

fprintf('Score was missing for %1.1f subjects \n',numel(toremove));
fprintf('R^2 is %1.2f with p value %1.2f \n',R_squared,p)
coefficients

%% B: Multiple linear regression: logistic regression
load('diagnoses.mat');
load('demographic.mat'); %added 2/13/2024
load('paper_indices.mat'); %added 2/13/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matter = 'WM'; %can change
VARIABLE = ARTICULATION; %ARTICULATION or STEREOTYPE (only n = 1, not enough for conclusions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(matter,'GM')
    load TBM_gray_down/PCA_TBM/Z96; 
end
if isequal(matter,'WM')
    load TBM_white_down/PCA_TBM/Z96; 
end
addpath 'Analysis codes'; 
[accuracy, kappa, sensitivity, specificity, PLDAprojected, w] = PLDA_knn(Z96, labels+2, 1, 3, 'brain');
close all;

% load the demographics
load('demographic'); 
load('paper_indices'); 
 
var = VARIABLE(inds==1); %articulation, coordination, enuresis, language, learning
projected_scores = PLDAprojected; 


%age, gender and volumes already proven to not be enough to differentiate
%the groups, not including in the model

%eliminate nans
toremove = find(isnan(var)); 
var(toremove) = []; projected_scores(toremove,:) = []; 

[B,dev,stats] = mnrfit(projected_scores,var+1,'Model','nominal'); 
[B stats.p]

%violin plots

L = labels; L(toremove) = []; %the labels corresponding to the included subjects

disorder = find(var==1); nodisorder = find(var==0);
Y{1} = projected_scores(disorder,1); 
Y{2} = projected_scores(nodisorder,1); 
violin(Y, 'xlabel',{'present','absent'},'facecolor',[0.5 0.5 0.5; 1 0.5 0.1],'facealpha',1); grid on; 
if isequal(matter,'WM')
    ylabel('WM TBM score'); 
elseif isequal(matter,'GM')
    ylabel('GM TBM score'); 
end
set(gca,'FontSize',17,'fontname','Times New Roman'); 
%set(gca,'fontname','Times New Roman');
%ylabel('TBM score')

%calculate sensitivity and specificity
C = confusionmat(var,sign(projected_scores(:,1)<0),'Order',[1 0]); %evaluate the degree to which negative projection score is sensitive for articulation disorder
sensitivity = C(1,1)./(C(1,1)+C(1,2))
specificity = C(2,2)./(C(2,2)+C(2,1))
%% What percentage of the articulation disorder subjects are deletion carriers vs. duplication carriers. 

load('demographic.mat')
load('diagnoses.mat')
load('paper_indices.mat')
indexes = find(inds==1); 

count = 0; N_del = 0; N_dup = 0; 
for i = 1:numel(indexes)
    if ARTICULATION(indexes(i))==1
        i
        if group(indexes(i))==-1
            N_del = N_del + 1;
        elseif group(indexes(i))==1
            N_dup = N_dup + 1; 
        end
       count = count + 1;
     end
end
N_dup./count
N_del./count



%% Convert .mat to .csv file
load('WM_TBM_visualization.mat');
test = image{3,5}; 
[X,Y,Z] = meshgrid(1:size(test,2),1:size(test,1),1:size(test,3));
A = [X(:) Y(:) Z(:) test(:)]; 


csvwrite('test_dup_1.csv', A,1,0);

%% Calculate volumes of the TBM-generated images for global tissue shifts analysis

thresh = 2.5; 
%gray matter only
load GM_TBM_visualization; 

dup = image{3,5}(36:end-35,26:end-25,26:end-25) > thresh*min(image{3,3}(36:end-35,26:end-25,26:end-25));
del = image{3,1}(36:end-35,26:end-25,26:end-25) > thresh*min(image{3,3}(36:end-35,26:end-25,26:end-25));  
con = image{3,3}(36:end-35,26:end-25,26:end-25) > thresh*min(image{3,3}(36:end-35,26:end-25,26:end-25));  
dup_half = image{3,4}(36:end-35,26:end-25,26:end-25) > thresh*min(image{3,3}(36:end-35,26:end-25,26:end-25));
del_half = image{3,2}(36:end-35,26:end-25,26:end-25) > thresh*min(image{3,3}(36:end-35,26:end-25,26:end-25));  

sum(dup(:))./numel(dup(:))
sum(dup_half(:))./numel(dup_half(:))
sum(con(:))./numel(con(:))
sum(del_half(:))./numel(del_half(:))
sum(del(:))./numel(del(:))

%white matter
load WM_TBM_visualization;
dup = image{3,5}(26:end-25,16:end-15,16:end-15) > thresh*min(image{3,3}(26:end-25,16:end-15,16:end-15));
del = image{3,1}(26:end-25,16:end-15,16:end-15) > thresh*min(image{3,3}(26:end-25,16:end-15,16:end-15));  
con = image{3,3}(26:end-25,16:end-15,16:end-15) > thresh*min(image{3,3}(26:end-25,16:end-15,16:end-15));  
dup_half = image{3,4}(26:end-25,16:end-15,16:end-15) > thresh*min(image{3,3}(26:end-25,16:end-15,16:end-15));
del_half = image{3,2}(26:end-25,16:end-15,16:end-15) > thresh*min(image{3,3}(26:end-25,16:end-15,16:end-15));  

sum(dup(:))./numel(dup(:))
sum(dup_half(:))./numel(dup_half(:))
sum(con(:))./numel(con(:))
sum(del_half(:))./numel(del_half(:))
sum(del(:))./numel(del(:))

%% %% Harvard-Oxford Atlas on displacement maps to determine areas of greatest shift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%matter = 'WM'; 
n = 3; %0 3 or -3 %standard deviations from the mean
WM_REG = [1 3 8 12 14]; %subcortical atlas values for white matter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:2
    if k==1
        matter = 'GM'; 
    elseif k==2
        matter = 'WM'; 
    end
%load the displacement fields
cd Visualization; 
load(strcat(matter,'_TBM_direction.mat')); 

%calculate the scores for standard deviations
sigma1 = std(PROJECTED(:,1)); 
sigma2 = std(PROJECTED(:,2)); 

%load the corresponding template image, in the same space as the template
if isequal(matter,'GM')
    load TBM_gray_down/T1_1_results.mat; 
    I0_trunc = imrotate(results.I0(26:end-25,36:end-35,26:end-25),90);
elseif isequal(matter,'WM')
    load TBM_white_down/T1_1_results.mat; 
    I0_trunc = imrotate(results.I0(16:end-15,26:end-25,16:end-15),90);
end
[M,N,K] =size(results.I0);
[X,Y,Z] = meshgrid(1:N,1:M,1:K); 

load jacobians2; %added just now
if k ==1
    load GM_TBM_visualization; 
    if n==-3 %added just now
        jacobian1 = GM_jacobian1_del;
        jacobian2 = GM_jacobian2_del; 
    elseif n==3
        jacobian1 = GM_jacobian1_dup;
        jacobian2 = GM_jacobian2_dup; 
    end

elseif k==2
    load WM_TBM_visualization; 

     if n==-3 %added just now
        jacobian1 = WM_jacobian1_del; 
        jacobian2 = WM_jacobian2_del; 
    elseif n==3
        jacobian1 = WM_jacobian1_dup; 
        jacobian2 = WM_jacobian2_dup; 
    end
end

%load atlases 
load atlases; %load atlas file
CORTICAL = round(cortical_atlas); SUBCORTICAL = round(subcortical_atlas); %round atlases to whole numbers
mask_cortical = (cortical_atlas~=0); %mask for gray matter images only 
mask_subcortical =(subcortical_atlas~=0); %mask for gray matter and white matter

%calculate normalized displacement maps
if isequal(matter,'WM')
    dmap1 = mask_subcortical.*jacobian1; 
    dmap2 = mask_subcortical.*jacobian2; 
 
    Zdmap1 = (dmap1-mean(dmap1(dmap1~=0)))./std(dmap1(dmap1~=0)); %zscore map
    Zdmap2 = (dmap2-mean(dmap2(dmap2~=0)))./std(dmap2(dmap2~=0)); %zscore map

    %subcortical atlas
    for i = 1:5
        mask = medfilt3(SUBCORTICAL==WM_REG(i)); %median filter to null edge signal
        SCD1(WM_REG(i)) = sum(Zdmap1(:).*mask(:))./nnz(mask); 
        SCD2(WM_REG(i)) = sum(Zdmap2(:).*mask(:))./nnz(mask);
        clear mask;
    end
elseif isequal(matter,'GM')
    dmap1 = mask_cortical.*jacobian1; 
    dmap2 = mask_cortical.*jacobian2; 
     
    Zdmap1 = (dmap1-mean(dmap1(dmap1~=0)))./std(dmap1(dmap1~=0)); %zscore map
    Zdmap2 = (dmap2-mean(dmap2(dmap2~=0)))./std(dmap2(dmap2~=0)); %zscore map

    %cortical atlas
    for i = 1:48
        mask = medfilt3(CORTICAL==i); %median filter to null edge signal
        CD1(i) = sum(sum(sum(Zdmap1.*mask)))./nnz(mask); %vector of normalized displacements for direction 1
        CD2(i) = sum(sum(sum(Zdmap2.*mask)))./nnz(mask); %vector of normalized displacmeents for direction 2
        clear mask; 
    end

    dmap1 = mask_subcortical.*jacobian1; 
    dmap2 = mask_subcortical.*jacobian2;

    Zdmap1 = (dmap1-mean(dmap1(dmap1~=0)))./std(dmap1(dmap1~=0)); %zscore map
    Zdmap2 = (dmap2-mean(dmap2(dmap2~=0)))./std(dmap2(dmap2~=0)); %zscore map
     %subcortical atlas
    for i = 1:21
        if ~any(WM_REG==i)
            mask = medfilt3(SUBCORTICAL==i); %median filter to null edge signal
            SCD1(i) = sum(sum(sum(Zdmap1.*mask)))./nnz(mask); %vector of normalized displacements for direction 1
            SCD2(i) = sum(sum(sum(Zdmap2.*mask)))./nnz(mask); %vector of normalized displacmeents for direction 2
            clear mask; 
        end
    end

end

end

round(CD1'*100)/100
[M,D] = sort(CD1,'descend'); 
[M' D']



