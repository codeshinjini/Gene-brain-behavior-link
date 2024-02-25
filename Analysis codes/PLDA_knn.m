%Shinjini Kundu (c) 2018
%Complete leave-three-out cross-validation 

function [accuracy, kappa, sensitivity, specificity, PLDAprojected, w] = PLDA_knn(X, labels, alpha, k, mode)
%This function performs classification in the discriminant subspace of the
%PLDA classifier 
%inputs:             
%                    k                number of tissue classes
%                    X                data matrix (n x d) 
%                    labels           labels of each sample
%                    alpha            penalization using the alpha parameter
%outputs:            accuracy         the accuracy of leave one out
%                                     cross-validation

fprintf('Newest version \n'); 
FINAL_FEATS = X; 
nPLDA = k-1; 
PERM_TESTS = 1000; %run 1000 permutation tests 
NUM_SUBJECTS = size(X,1); 
%NEIGHBORS = round(sqrt(NUM_SUBJECTS)); %sqrt of the number of data points

%visualize how far apart the clusters lie
figure;
if(isempty(alpha))
    fprintf('Now computing alpha \n'); 
    Curveoption.low = 0; %lower bound of alpha 0.01
    Curveoption.high = 0.5; %upper bound of alpha
    Curveoption.step = 0.1; %step size of increaseing alpha
    Curveoption.nPLDA = nPLDA; %the dimension of subspaces to be compared
    alpha = Calculate_Alpha( FINAL_FEATS, labels, Curveoption, 1 );
    [Vec, eig] = PLDA(FINAL_FEATS', labels, alpha, nPLDA);
else
    [Vec, eig] = PLDA(FINAL_FEATS', labels, alpha, nPLDA);
end

% Make sure that sign of direction is consistent, for deletions(-) and
% control (+)
PLDAprojected = FINAL_FEATS*Vec;
if sum(sign(PLDAprojected((labels==1),1)))>0 && sum(sign(PLDAprojected(labels==3,1)))<0
    Vec(:,1) = -Vec(:,1); 
    fprintf('switch direction 1\n')
end
if sum(sign(PLDAprojected(labels==2,2)))<0
    Vec(:,2) = -Vec(:,2);
    fprintf('switch direction 2\n')
end

%After standardization, compute PLDAprojected
PLDAprojected = FINAL_FEATS*Vec;

% % % orange = [1 0.5 0.2];
% % % h = plot(PLDAprojected(labels==-1,1),PLDAprojected(labels==-1,2),'.k','MarkerSize',25);
% % % set(h,'Color',orange);
% % % hold on;
% % % plot(PLDAprojected(labels==1,1),PLDAprojected(labels==1,2),'.b','MarkerSize',25); 
% % % h2 = plot(PLDAprojected(labels==0,1),PLDAprojected(labels==0,2),'.g','MarkerSize',25); 
% % % %gree = [0 1 0]; 
% % % %set(h2,'Color',gree); 
% % % set(gca,'FontSize',18); 
% % % set(gca,'FontName','Times New Roman'); 
% % % xlabel('Discriminant direction 1'); 
% % % ylabel('Discriminant direction 2'); 

color(1,:) = [0.5 0.5 0.5]; color(2,:) = [0 0.4 0.4]; color(3,:) = [1 0.5 0.1]; 

for i = 1:k
   h = plot(PLDAprojected(labels==i,1),PLDAprojected(labels==i,2),'.k','MarkerSize',25);
   set(h,'Color',color(i,:));
   hold all;
end
%make the axes be in terms of standard deviations of projection score for
%each axis rather than numbers themselves
mean1 = mean(PLDAprojected(:,1)); sigma1 = std(PLDAprojected(:,1)); 
mean2 = mean(PLDAprojected(:,2)); sigma2 = std(PLDAprojected(:,2)); 
xl = xlim; yl = ylim;
if isequal(mode,'brain')
    xticks([xl(1) -1.5*sigma1 0 1.5*sigma1 xl(2)]+mean1); 
    xticklabels({'','-1.5\sigma_1','0','1.5\sigma_1',''});
    try
        yticks([yl(1) -2*sigma2 -sigma2 0 sigma2 2*sigma2 yl(2)]+mean2); 
        yticklabels({'','-2\sigma_2','','0','','2\sigma_2',''});
    catch 
        yticks([yl(1) -1.5*sigma2 -sigma2 0 sigma2 1.5*sigma2 yl(2)]+mean2); 
        yticklabels({'','-1.5\sigma_2','','0','','1.5\sigma_2',''});
    end
else
    xticks([xl(1) mean1-sigma1 mean1 mean1+sigma1 xl(2)]); 
    xticklabels({'','-\sigma_1','0','\sigma_1',''});
    yticks([yl(1) mean2-sigma2 mean2 mean2+sigma2  yl(2)]); 
    yticklabels({'','-\sigma_2','0','\sigma_2',''});
end

legend('deletion','control','duplication'); 
set(gca,'FontSize',18); 
set(gca,'FontName','Times New Roman'); 
xlabel('Discriminant direction 1'); 
ylabel('Discriminant direction 2'); 
grid on;
%%%%%%%%%%%%%

%Perform centroid classification using leave-one-out cross-validation 
fprintf('Now computing classifier and boundary \n'); 
%find the number of unique vectors
% classes = unique(labels);
% for i = 1:numel(classes)
%     classInds{i} = find(labels==classes(i)); 
% end
% combinations = combvec(classInds{1},classInds{2}, classInds{3}); 
 
w = 0; %weighted direction 
for i = 1:size(X,1)
  i
  test_labels = labels(i); 
  X_train = FINAL_FEATS; train_labels = labels; X_train(i,:) = []; train_labels(i) = []; 
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
 
  %make sure the directions are standardized with control (+) and deletions
  %(-1)
  pr = X_train*Vec_train; 
  if sum(sign(pr(train_labels==1,1))) >0  && sum(sign(pr(train_labels==3,1)))<0
    Vec_train(:,1) = -Vec_train(:,1); 
  end
  if sum(sign(pr(train_labels==2,2))) <0 
    Vec_train(:,2) = -Vec_train(:,2);
  end
    
  %now compute the projections
  X_test = FINAL_FEATS(i,:); x_tilde = X_test*Vec_train;
  w = w + Vec_train; 
  
  %assign labels to test data
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
  predicted(i,:) = p;
  actual(i,:) = test_labels';
end

%compute the information needed for the results

w = w./size(X,1); %average of all classification boundaries

%compute overall test accuracy and cohen's kappa
t = (actual(:)~=predicted(:)); 
cms{1} = confusionmat(actual(:),predicted(:));
kappa = CohensKappa(cms); 
accuracy  = 1-sum(t)/numel(actual(:));

%calculating sensitivity and specificity 
for i = 1:k
    truth = actual(:); truth = (truth==i);  
    pred = predicted(:); pred =(pred==i);
    CP = classperf(truth,pred); 
    sensitivity(i) = CP.Sensitivity; 
    specificity(i) = CP.Specificity; 
    %[C, order] = confusionmat(truth, pred); 
    %posInd = find(order == 1); 
    %negInd = find(order == 0); 
    %sensitivity(i) = C(posInd, posInd)/(C(posInd,posInd)+C(posInd,negInd)); 
    %specificity(i) = C(negInd,negInd)/(C(negInd,negInd)+C(negInd,posInd));       
end

boundary = 0;
clear distances; 
%visualize the classifier boundaries
[Xv,Yv] = meshgrid(xl(1):10:xl(2),yl(1):10:yl(2));
for j = 1:size(Xv,1)
    for l = 1:size(Xv,2)
        for m = 1:k
            distances(m,j,l) = norm([Xv(j,l) Yv(j,l)]-mu_hat(m,:),2)^2-log(pi(m));
        end
        [~,p] = min(distances(:,j,l)); 
        classes(j,l) = p; 
    end
end
boundary = boundary + classes;
figure; 
%markers = ['*', 's','o']; 
for i = 1:k
   %h = plot(PLDAprojected(labels==i,1),PLDAprojected(labels==i,2),markers(i),'linewidth',2);
   h = plot(PLDAprojected(labels==i,1),PLDAprojected(labels==i,2),'.k','MarkerSize',25);
   set(h,'Color',color(i,:));
   hold all;
end
%make the axes be in terms of standard deviations of projection score for
%each axis rather than numbers themselves
mean1 = mean(PLDAprojected(:,1)); sigma1 = std(PLDAprojected(:,1)); 
mean2 = mean(PLDAprojected(:,2)); sigma2 = std(PLDAprojected(:,2)); 
if isequal(mode,'brain')
    xticks([xl(1) -1.5*sigma1 0 1.5*sigma1 xl(2)]+mean1); 
    xticklabels({'','-1.5\sigma_1','0','1.5\sigma_1',''});
    try
        yticks([yl(1) -2*sigma2 -sigma2 0 sigma2 2*sigma2 yl(2)]+mean2); 
        yticklabels({'','-2\sigma_2','','0','','2\sigma_2',''});
    catch 
        yticks([yl(1) -1.5*sigma2 -sigma2 0 sigma2 1.5*sigma2 yl(2)]+mean2); 
        yticklabels({'','-1.5\sigma_2','','0','','1.5\sigma_2',''});
    end
else
    xticks([xl(1) mean1-sigma1 mean1 mean1+sigma1 xl(2)]); 
    xticklabels({'','-\sigma_1','0','\sigma_1',''});
    yticks([yl(1) mean2-sigma2 mean2 mean2+sigma2  yl(2)]); 
    yticklabels({'','-\sigma_2','0','\sigma_2',''});
end

set(gca,'FontSize',18); 
set(gca,'FontName','Times New Roman'); 
xlabel('Discriminant direction 1'); 
ylabel('Discriminant direction 2'); 
grid on;
contour(Xv,Yv,round(boundary)-1,[1 2],'k','LineWidth',2); %,[0.5 0.5],'k'); 
legend('deletion','control','duplication'); 

% % % % % fprintf('Now running permutation tests \n')
% % % % % %compute significance level using permutation tests
% % % % % if (permtest ==1)
% % % % %     p_value = 0; 
% % % % %     for m = 1:PERM_TESTS
% % % % %         fprintf('Now on the %d th test \n', m); 
% % % % %         permuted = labels(randperm(NUM_SUBJECTS)); 
% % % % %         for i = 1:size(combinations,2)
% % % % %         test_labels = permuted(combinations(:,i)); 
% % % % %         X_train = FINAL_FEATS; train_labels = permuted; X_train(combinations(:,i),:) = []; train_labels(combinations(:,i)) = []; 
% % % % %         %Calcualte alpha
% % % % %         if(isnan(alpha))
% % % % %             Curveoption.low = 0; %lower bound of alpha 0.01
% % % % %             Curveoption.high = 0.5; %upper bound of alpha
% % % % %             Curveoption.step = 0.1; %step size of increaseing alpha
% % % % %             Curveoption.nPLDA = nPLDA; %the dimension of subspaces to be compared
% % % % %             alpha = Calculate_Alpha( X_train, train_labels, Curveoption,0 );
% % % % %             [Vec_train, ~] = PLDA(X_train', train_labels, alpha, nPLDA);
% % % % %             alpha = nan;
% % % % %         else
% % % % %             [Vec_train, ~] = PLDA(X_train', train_labels, alpha, nPLDA);
% % % % %         end
% % % % %  
% % % % %         X_test = FINAL_FEATS(combinations(:,i),:); 
% % % % %         x_tilde = X_test*Vec_train;
% % % % %   
% % % % %         %assign labels to test data
% % % % %         for j= 1:k
% % % % %             pi(j) = sum(train_labels==j)/numel(train_labels);  
% % % % %             mu(j,:) = sum(X_train((train_labels==j),:))./sum(train_labels==j);
% % % % %             mu_hat(j,:) = mu(j,:)*Vec_train;
% % % % %             for l = 1:size(x_tilde,1)
% % % % %                 distances(:,j,l) = norm(x_tilde(l,:)-mu_hat(j,:),2)-log(pi(j));
% % % % %             end
% % % % %         end
% % % % %         for l = 1:size(x_tilde,1)
% % % % %             [~,p(l)] = min(distances(:,:,l));
% % % % %         end 
% % % % %         predicted(i,:) = p;
% % % % %         actual(i,:) = test_labels';
% % % % %         end
% % % % %     end
% % % % %     ptest_accuracy  = 1-sum(t)/numel(actual(:));
% % % % %     if ptest_accuracy > accuracy
% % % % %         p_value = p_value + 1;
% % % % %     end
% % % % %     p_value = p_value/PERM_TEST; 
% % % % % else
% % % % %     p_value = nan; 
% % % % % end

end

