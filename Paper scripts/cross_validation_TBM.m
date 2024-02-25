%number of cross validation iterations

function cross_validation_TBM(matter)

T = 1000; %number of cross-validation iterations 

%cluster

if isequal(matter,'WM')
    cd classifications/WM/
elseif isequal(matter,'GM')
    cd classifications/GM/
elseif isequal(matter,'WM_unmodulated')
    cd classifications/WM_unmodulated/
elseif isequal(matter,'GM_unmodulated')
    cd classifications/GM_unmodulated/
elseif isequal(matter,'WM_modulated')
    cd classifications/WM_modulated/
elseif isequal(matter,'GM_modulated')
    cd classifications/GM_modulated/
elseif isequal(matter,'DBM')
    cd classifications/DBM/
elseif isequal(matter,'WM_orig')
    cd classifications/WM_orig/
elseif isequal(matter,'GM_orig')
    cd classifications/GM_orig/
elseif isequal(matter,'cov_WM')
    cd classifications/cov_WM_TBM/
elseif isequal(matter,'cov_GM')
    cd classifications/cov_GM_TBM/
elseif isequal(matter,'Tot')
    cd classifications/Tot/
end

for j = 1:T
    load(strcat('iter',num2str(j))); 
    ACC(j) = accuracy; KAPPA(j) = kappa; 
    SENS(j,:) = sensitivity; SPEC(j,:) = specificity; 
%    ALPHAS(j) = alpha; 
    GROUND_TRUTH(j,:) = actual; 
    PREDICTED(j,:) = predicted; 
    clear accuracy kappa sensitivity specificity actual predicted; 
    j
end

[average, lowerCI,higherCI] = classification_compute_CI(ACC,T); 
fprintf('The mean accuracy is %0.1f with CI between %0.1f and %0.1f \n',average*10^2, lowerCI*10^2,higherCI*10^2); 

[average, lowerCI,higherCI] = classification_compute_CI(KAPPA,T); 
fprintf('The mean kappa is %0.2f with CI between %0.2f and %0.2f \n',average, lowerCI,higherCI); 

for k = 1:3
    [average, lowerCI,higherCI] = classification_compute_CI(SENS(:,k),T); 
    fprintf('The mean sensitivity is %0.1f with CI between %0.1f and %0.1f \n',average*10^2, lowerCI*10^2,higherCI*10^2); 
end

for k = 1:3
    [average, lowerCI,higherCI] = classification_compute_CI(SPEC(:,k),T); 
    fprintf('The mean specificity is %0.1f with CI between %0.1f and %0.1f \n',average*10^2, lowerCI*10^2,higherCI*10^2); 
end


end