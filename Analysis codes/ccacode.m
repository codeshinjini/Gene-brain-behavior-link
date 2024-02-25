load paper_indices; load demographic; 
LABELS = group(inds==1);
del_inds = find(LABELS==-1); dup_inds = find(LABELS==1); con_inds = find(LABELS==0); 

%load the feature matrices
load TBM_gray_down/features;
GM = features; 
load TBM_white_down/features;
WM = features;

%%things to modify
GRP_INDS = 1:206;
numIter = 100; 
%%%%%%%%%%%%%%%%%%

N_train = round(0.75*numel(GRP_INDS));
N_test = numel(GRP_INDS)-N_train;

r_test = 0; r_train = 0; 
averageA = 0; averageB = 0; 

rng(1); %for reproducibility
for i = 1:numIter
    i
    indices = randsample(numel(GRP_INDS),N_train); 
    X_train = GM(GRP_INDS(indices),:); Y_train = WM(GRP_INDS(indices),:); 
    X_test = GM(GRP_INDS,:); X_test(indices,:) = []; Y_test = WM(GRP_INDS,:); Y_test(indices,:) = [];  

    %run PCA
    [~,~,~,EIGENV_X,Z_train_X,Z_test_X] = PCA_decomp(X_train, X_test);
    [~,~,~,EIGENV_Y,Z_train_Y,Z_test_Y] = PCA_decomp(Y_train, Y_test);

    %run CCA
    [A,B,r,U_train,V_train,stats] = canoncorr(Z_train_X,Z_train_Y) ;
    U_test = (Z_test_X)*A; 
    V_test = (Z_test_Y)*B;
    R = corrcoef(U_test(:,1),V_test(:,1)); %first canonical variable

    %project direction back onto TBM space
    %only look at the first CCA direction
    
    averageA = EIGENV_X*A(:,1) + averageA; 
    averageB = EIGENV_Y*B(:,1) + averageB; 

   %calculate r_train and r_test
    r_train = r(1) + r_train; 
    r_test = R(2,1) + r_test; 
    SCATTER_UV{i} = [U_test(:,1) V_test(:,1)]; 
end
r_test = r_test/numIter; %test correlation coefficient
averageA = averageA/numIter; averageB = averageB/numIter; 
r_train = r_train/numIter; 

% visualizing the directions found by CCA
 
 %U_vis = (GM-repmat(mean(GM),206,1))*averageA; 
 %V_vis = (WM-repmat(mean(WM),206,1))*averageB;
 %R = corrcoef(U_vis(:,1),V_vis(:,1)); 
save CCA_results r_test r_train SCATTER_UV

%% scatter plot
%%%%%%%%%%%%%can change%%%%%%%%%%%%%%%%%%%
numIter = 100; p = 0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_vis = []; V_vis = []; 
for i = 1:numIter
    U_vis = [U_vis;SCATTER_UV{i}(:,1)]; 
    V_vis = [V_vis;SCATTER_UV{i}(:,2)]; 
end
plot(U_vis(:,1),V_vis(:,1),'.','MarkerSize',3); grid on; 
xlabel('projection of GM'); ylabel('projection of WM');  
sigma_U = std(U_vis(:,1)); sigma_V = std(V_vis(:,1));  
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',18);
set(gca,'YTick',[-2*sigma_V, -sigma_V, 0, sigma_V, 2*sigma_V]);
set(gca,'YTickLabel',[{'-2\sigma_2','-\sigma_2','0','\sigma_2','2\sigma_2'}]);
set(gca,'XTick',[-2*sigma_U, -sigma_U, 0, sigma_U, 2*sigma_U]);
set(gca,'XTickLabel',[{'-2\sigma_1','-\sigma_1','0','\sigma_1','2\sigma_1'}]); 

dim = [.2 .5 .3 .3];
if p==0
    str = {strcat('R = ',num2str(round(r_test*10^2)/10^2)), strcat('p <0.01 ')}; %= ',num2str(results.p))};
else
    str = {strcat('R = ',num2str(round(results.r_test*10^2)/10^2)), strcat('p = ',num2str(round(p*10^2)/10^2))}; %= ',num2str(results.p))};
end
annotation('textbox',dim,'String',str,'position',[0.73 0.25,0.1,0.1],'FitBoxToText','on');

%% permutation testing
r_test_perm = 0; r_train_perm = 0; p = 0;
for j = 1:100 %number of permutation tests
    permuted = randsample(numel(GRP_INDS),numel(GRP_INDS));
    WM_permuted = WM(permuted,:);
    j
    rng(1);
    for i = 1:numIter
        i
    	indices = randsample(numel(GRP_INDS),N_train); 
    	X_train = GM(GRP_INDS(indices),:); Y_train = WM_permuted(GRP_INDS(indices),:); 
    	X_test = GM(GRP_INDS,:); X_test(indices,:) = []; Y_test = WM_permuted(GRP_INDS,:); Y_test(indices,:) = [];  

    	%run PCA
    	[~,~,~,EIGENV_X,Z_train_X,Z_test_X] = PCA_decomp(X_train, X_test);
    	[~,~,~,EIGENV_Y,Z_train_Y,Z_test_Y] = PCA_decomp(Y_train, Y_test);

    	%run CCA
    	[A,B,r,U_train,V_train,stats] = canoncorr(Z_train_X,Z_train_Y) ;
    	U_test = (Z_test_X)*A; 
    	V_test = (Z_test_Y)*B;
    	R = corrcoef(U_test(:,1),V_test(:,1)); %first canonical variable

    	%project direction back onto TBM space
   	 %only look at the first CCA direction
    
   	 averageA = EIGENV_X*A(:,1) + averageA; 
   	 averageB = EIGENV_Y*B(:,1) + averageB; 

   	%calculate r_train and r_test
    	r_train_perm = r(1) + r_train_perm; 
    	r_test_perm = R(2,1) + r_test_perm; 

    end
    r_test_perm = r_test_perm/numIter; 
    r_train_perm = r_train_perm/numIter;
    if r_test_perm >= r_test
        p = p + 1;
    end
end
p_value = p/100
 