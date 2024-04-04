%Shinjini Kundu (c) 2018
%Complete cross-validation 

function [accuracy, kappa, sensitivity, specificity, w,alpha,actual,predicted] = PLDA_knn_CV(X, labels, alpha, k, mode,ITER)
%This function performs classification in the discriminant subspace of the
%PLDA classifier 
%inputs:             
%                    k                number of tissue classes
%                    X                data matrix (n x d) %full data matrix
%                    before PCA
%                    labels           labels of each sample
%                    alpha            penalization using the alpha parameter
%outputs:            accuracy         the accuracy of leave one out
%                                     cross-validation

VIZ = 0; %change to 1 if calculate visualization 
load CVpartitions; %load the cross-validation partitions
FINAL_FEATS = X; 
c = C{ITER}; %load the specific partition of interest
folds = 10; 

nPLDA = 2; 
%NUM_SUBJECTS = size(X,1); 

%sc = 0; %to display scatter plot or not
% % % % 
% % % % %visualize how far apart the clusters lie
% % % % figure;
% % % % if(isempty(alpha))
% % % %     fprintf('Now computing alpha \n'); 
% % % %     Curveoption.low = 0; %lower bound of alpha 0.01
% % % %     Curveoption.high = 0.5; %upper bound of alpha
% % % %     Curveoption.step = 0.1; %step size of increaseing alpha
% % % %     Curveoption.nPLDA = nPLDA; %the dimension of subspaces to be compared
% % % %     alpha = Calculate_Alpha( FINAL_FEATS, labels, Curveoption, 1 );
% % % %     [Vec, eig] = PLDA(FINAL_FEATS', labels, alpha, nPLDA);
% % % % else
% % % %     [Vec, eig] = PLDA(FINAL_FEATS', labels, alpha, nPLDA);
% % % % end
% % % % PLDAprojected = FINAL_FEATS*Vec;

% % % % % % orange = [1 0.5 0.2];
% % % % % % h = plot(PLDAprojected(labels==-1,1),PLDAprojected(labels==-1,2),'.k','MarkerSize',25);
% % % % % % set(h,'Color',orange);
% % % % % % hold on;
% % % % % % plot(PLDAprojected(labels==1,1),PLDAprojected(labels==1,2),'.b','MarkerSize',25); 
% % % % % % h2 = plot(PLDAprojected(labels==0,1),PLDAprojected(labels==0,2),'.g','MarkerSize',25); 
% % % % % % %gree = [0 1 0]; 
% % % % % % %set(h2,'Color',gree); 
% % % % % % set(gca,'FontSize',18); 
% % % % % % set(gca,'FontName','Times New Roman'); 
% % % % % % xlabel('Discriminant direction 1'); 
% % % % % % ylabel('Discriminant direction 2'); 
% % % 
% % % color(1,:) = [0.5 0.5 0.5]; color(2,:) = [0 0.4 0.4]; color(3,:) = [1 0.5 0.1]; 
% % % 
% % % for i = 1:k
% % %    h = plot(PLDAprojected(labels==i,1),PLDAprojected(labels==i,2),'.k','MarkerSize',25);
% % %    set(h,'Color',color(i,:));
% % %    hold all;
% % % end
% % % %make the axes be in terms of standard deviations of projection score for
% % % %each axis rather than numbers themselves
% % % mean1 = mean(PLDAprojected(:,1)); sigma1 = std(PLDAprojected(:,1)); 
% % % mean2 = mean(PLDAprojected(:,2)); sigma2 = std(PLDAprojected(:,2)); 
% % % xl = xlim; yl = ylim;
% % % if isequal(mode,'brain')
% % %     xticks([xl(1) -1.5*sigma1 0 1.5*sigma1 xl(2)]+mean1); 
% % %     xticklabels({'','-1.5\sigma_1','0','1.5\sigma_1',''});
% % %     yticks([yl(1) -2*sigma2 -sigma2 0 sigma2 2*sigma2 yl(2)]+mean2); 
% % %     yticklabels({'','-2\sigma_2','','0','','2\sigma_2',''});
% % % else
% % %     xticks([xl(1) mean1-sigma1 mean1 mean1+sigma1 xl(2)]); 
% % %     xticklabels({'','-\sigma_1','0','\sigma_1',''});
% % %     yticks([yl(1) mean2-sigma2 mean2 mean2+sigma2  yl(2)]); 
% % %     yticklabels({'','-\sigma_2','0','\sigma_2',''});
% % % end
% % % 
% % % legend('deletion','control','duplication'); 
% % % set(gca,'FontSize',18); 
% % % set(gca,'FontName','Times New Roman'); 
% % % xlabel('Discriminant direction 1'); 
% % % ylabel('Discriminant direction 2'); 
% % % grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Perform centroid classification using leave-one-out cross-validation 
%fprintf('Now computing classifier and boundary \n'); 
%find the number of unique vectors
% % classes = unique(labels);
% % for i = 1:numel(classes)
% %     classInds{i} = find(labels==classes(i)); 
% % end
%% combinations = [1:205]; %indices that depend on cross validation procedure
 
w = []; %direction variable
dir1 = 0; dir2 = 0; %weighted directions
MEAN = 0; %average of A_means
predicted = []; actual = [];
for i = 1:folds
  fprintf('Running fold %d \n',i);
  test_labels = labels(test(c,i));
  if size(FINAL_FEATS,2) < size(FINAL_FEATS,1)
    X_train = FINAL_FEATS; train_labels = labels; X_train(test(c,i),:) = []; train_labels(test(c,i)) = [];
  elseif size(FINAL_FEATS,2) > size(FINAL_FEATS,1)
    train_features = FINAL_FEATS(training(c,i),:); train_labels = labels(training(c,i));
    test_features = FINAL_FEATS(test(c,i),:); test_labels = labels(test(c,i));  
    [Z96_train,Z96_test,A_mean,EIGENV] = PCA_decomp(train_features, test_features); 
    X_train = Z96_train; 
  end
  %Calcualte alpha
  if(isnan(alpha))
    Curveoption.low = 0; %lower bound of alpha 0.01
    Curveoption.high = 0.5; %upper bound of alpha
    Curveoption.step = 0.1; %step size of increaseing alpha
    Curveoption.nPLDA = nPLDA; %the dimension of subspaces to be compared
    alpha = Calculate_Alpha( X_train, train_labels, Curveoption,0 );
    [Vec_train, ~] = PLDA(X_train', train_labels, alpha, nPLDA);
    alpha = nan;
  else
    [Vec_train, ~] = PLDA(X_train', train_labels, alpha, nPLDA);
  end
 
  if size(FINAL_FEATS,2) < size(FINAL_FEATS,1)
      X_test = FINAL_FEATS(test(c,i),:); 
  elseif size(FINAL_FEATS,2) > size(FINAL_FEATS,1)
      X_test =  Z96_test; 
  end
  x_tilde = X_test*Vec_train; %PLDAprojected

  %% for visualization
  if VIZ
    v1 = Vec_train(:,1); %x direction 
    v2 = Vec_train(:,2); %y direction

    V1 = zeros(size(X,2),1); V2 = V1; 
    for m = 1:size(X_train,2)
        V1 = V1 + v1(m)*EIGENV(:,m); 
        V2 = V2 + v2(m)*EIGENV(:,m);
    end
    % add the directions V1 and V2 to the running total across each fold 
    dir1 = dir1 + V1;
    dir2 = dir2 + V2; 
    MEAN = MEAN + A_mean; 
  end

  %assign labels to test data
  clear pi mu mu_hat
  for j = 1:k
    pi(j) = sum(train_labels==j)/numel(train_labels);  
    mu(j,:) = sum(X_train((train_labels==j),:))./sum(train_labels==j);
    mu_hat(j,:) = mu(j,:)*Vec_train;
    for l = 1:size(x_tilde,1)
        distances(j,l) = 0.5*norm(x_tilde(l,:)-mu_hat(j,:),2)^2-log(pi(j));
    end
  end
  for l = 1:size(x_tilde,1)
    [~,p(l)] = min(distances(:,l));
  end 
  predicted = [predicted p];
  actual = [actual test_labels];
  clear distances p; 
end

%compute the information needed for the results

if VIZ
    w.dir1 = dir1./folds; %average of all classification boundaries
    w.dir2 = dir2./folds; 
    w.MEAN = MEAN./folds; 
end

%compute overall test accuracy and cohen's kappa
t = (actual(:)~=predicted(:)); 
cms{1} = confusionmat(actual(:),predicted(:));
kappa = CohensKappa(cms); 
accuracy  = 1-sum(t)/numel(actual(:));

%calculating sensitivity and specificity 
for i = 1:k
    truth = actual(:); truth = (truth==i);  
    pred = predicted(:); pred =(pred==i);
    %CP = classperf(truth,pred); 
    %sensitivity(i) = CP.Sensitivity; 
    %specificity(i) = CP.Specificity; 
    [C, order] = confusionmat(truth, pred); 
    posInd = find(order == 1); 
    negInd = find(order == 0); 
    sensitivity(i) = C(posInd, posInd)/(C(posInd,posInd)+C(posInd,negInd)); 
    specificity(i) = C(negInd,negInd)/(C(negInd,negInd)+C(negInd,posInd));       
end

% % if sc ==1
% % boundary = 0;
% % clear distances; 
% % %visualize the classifier boundaries
% % [Xv,Yv] = meshgrid(xl(1):10:xl(2),yl(1):10:yl(2));
% % for j = 1:size(Xv,1)
% %     for l = 1:size(Xv,2)
% %         for m = 1:k
% %             distances(m,j,l) = norm([Xv(j,l) Yv(j,l)]-mu_hat(m,:),2)^2-log(pi(m));
% %         end
% %         [~,p] = min(distances(:,j,l)); 
% %         classes(j,l) = p; 
% %     end
% % end
% % boundary = boundary + classes;
% % figure; 
% % for i = 1:k
% %    h = plot(PLDAprojected(labels==i,1),PLDAprojected(labels==i,2),'.k','MarkerSize',25);
% %    set(h,'Color',color(i,:));
% %    hold all;
% % end
% % %make the axes be in terms of standard deviations of projection score for
% % %each axis rather than numbers themselves
% % mean1 = mean(PLDAprojected(:,1)); sigma1 = std(PLDAprojected(:,1)); 
% % mean2 = mean(PLDAprojected(:,2)); sigma2 = std(PLDAprojected(:,2)); 
% % if isequal(mode,'brain')
% %     xticks([xl(1) -1.5*sigma1 0 1.5*sigma1 xl(2)]+mean1); 
% %     xticklabels({'','-1.5\sigma_1','0','1.5\sigma_1',''});
% %     yticks([yl(1) -2*sigma2 -sigma2 0 sigma2 2*sigma2 yl(2)]+mean2); 
% %     yticklabels({'','-2\sigma_2','','0','','2\sigma_2',''});
% % else
% %     xticks([xl(1) mean1-sigma1 mean1 mean1+sigma1 xl(2)]); 
% %     xticklabels({'','-\sigma_1','0','\sigma_1',''});
% %     yticks([yl(1) mean2-sigma2 mean2 mean2+sigma2  yl(2)]); 
% %     yticklabels({'','-\sigma_2','0','\sigma_2',''});
% % end
% % 
% % set(gca,'FontSize',18); 
% % set(gca,'FontName','Times New Roman'); 
% % xlabel('Discriminant direction 1'); 
% % ylabel('Discriminant direction 2'); 
% % grid on;
% % contour(Xv,Yv,round(boundary)-1,[1 2],'k','LineWidth',2); %,[0.5 0.5],'k'); 
% % legend('deletion','control','duplication');
% % end
% % 

end

