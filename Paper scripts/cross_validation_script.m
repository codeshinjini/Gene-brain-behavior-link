%number of cross validation iterations

function cross_validation_script(features,label)

T = 1000; %number of cross-validation iterations 

%cluster

for j = 1:T
    [accuracy, kappa, sensitivity, specificity, w, alpha, actual, predicted] = PLDA_knn_CV(features, label+2, 1, 3,'features', j); %run this on the
    ACC(j) = accuracy; KAPPA(j) = kappa; 
    SENS(j,:) = sensitivity; SPEC(j,:) = specificity; 
    GROUND_TRUTH(j,:) = actual; 
    PREDICTED(j,:) = predicted; 
    clear accuracy kappa sensitivity specificity actual predicted; 
    j
end

[average, lowerCI,higherCI] = classification_compute_CI(ACC,T); 
fprintf('The mean accuracy is %f with CI between %f and %f \n',average*10^2, lowerCI*10^2,higherCI*10^2); 

[average, lowerCI,higherCI] = classification_compute_CI(KAPPA,T); 
fprintf('The mean kappa is %f with CI between %f and %f \n',average, lowerCI,higherCI); 

for k = 1:3
    [average, lowerCI,higherCI] = classification_compute_CI(SENS(:,k),T); 
    fprintf('The mean sensitivity is %f with CI between %f and %f \n',average*10^2, lowerCI*10^2,higherCI*10^2); 
end

for k = 1:3
    [average, lowerCI,higherCI] = classification_compute_CI(SPEC(:,k),T); 
    fprintf('The mean specificity is %f with CI between %f and %f \n',average*10^2, lowerCI*10^2,higherCI*10^2); 
end



end