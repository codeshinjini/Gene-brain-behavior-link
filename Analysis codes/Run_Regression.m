%Shinjini Kundu (c) 2016
%Note: when using covariates, check to make sure that scaling of each
%column is in the same order of magnitude as the other columns. If not,
%must whiten the covariates matrix (column mean 0 and standard deviation 1)

function [ results ] = Run_Regression( matrix, y, frac, I0, A_mean, mode, visualize, covariates )
%Finds regression direction in the TBM space or image space that correlates 
%most with the independent variable of interest. Code also can adjust to 
%remove the effects of confounding variables
%inputs:   matrix              Data matrix, oriented NxD, where N is
%                               subjects (assumes that this is PCA-reduced
%                               data matrix). For original data matrix,
%                               must modify code to mean-subtract so means
%                               of each column are 0
%          y                   independent variable in column matrix
%          frac                number between 0 and 1 indicating fraction of data for training
%          I0                  smoothed template from which TBM images can be generated
%          A_mean              mean of feature matrix in PCA
%          mode                'TBM' or 'image'
%          visualize           0 or 1, depending on if images are desired to be generated
%          covariates          matrix of confounding variables, oriented
%                               Nxd, where 
%output:   wcorr               most correlated direction in the space
%          image               cell matrix of images corresponding to n
%          lambda              std of projection scores
%          CC_train            Pearson's correlation coefficient for
%                               training set
%          CC_test             Pearson's correlation coefficient for
%                               testing set
%          scores              projection scores
%          train_inds          indices used for training
%          test_inds           indices used for testing
%          p                   p value computed by permutation testing
%          n                   n std


[M,N,K] = size(I0); 
results.I0 = I0; results.A_mean = A_mean;

[X,Y,Z] = meshgrid(1:N,1:M,1:K); 
NUM_TESTS = 10000; 
close all; 
Data = matrix'; %orient matrix DxN
%X_mean = mean(Data')'; %don't need to subtract mean for the PCA space
num = size(Data,2);
vizaxis = 'x'; %determine whether you want to visualize brains acros the projection score(y) or independent variable(x)
%Data = Data-repmat(X_mean,1,num); %mean-subtracted matrix with images as columns

% % % %cross-validation code
% % % for i = 1:num
% % %     test_inds = i; 
% % %     train_inds = 1:num; 
% % %     train_inds(test_inds) = []; 
    
    %create testing and training indices    
    train_inds = randsample(num,round(num*frac),'false'); %randomly select 70% of samples for training
    test_inds = 1:num; 
    test_inds(train_inds) = []; 
        
    if ~isempty(covariates)
        v = y - covariates*((covariates'*covariates)\covariates')*y; 
    elseif isempty(covariates)
        v = y;
    end

    
    %%%%%Debugging code to make sure residual v and LS are actually orthogonal
%     LS = covariates/(covariates'*covariates)*covariates'*y;
%     CosTheta = max(min(dot(v,LS)/(norm(v)*norm(LS)),1),-1);
% ThetaInDegrees = real(acosd(CosTheta))
%%%%%
    %wcorr = Data(:,train_inds)*v(train_inds)/sqrt(v(train_inds)'*Data(:,train_inds)'*Data(:,train_inds)*v(train_inds)); %most correlated direction 
    wcorr = Data(:,train_inds)*(v(train_inds)-mean(v(train_inds)))/sqrt((v(train_inds)-mean(v(train_inds)))'*Data(:,train_inds)'*Data(:,train_inds)*(v(train_inds)-mean(v(train_inds))));
    lambda = std(Data'*wcorr); 


    Xw_train=Data(:,train_inds)'*wcorr; %Projection of the dataset onto the most correlative direction 
    Xw_test=Data(:,test_inds)'*wcorr;
    
    %scores(i) = Xw_test; 
    scores(train_inds) = Xw_train; scores(test_inds) = Xw_test; 
    
   
    %Calculate the correlation coefficient
    
        CC_train = corr(double(Xw_train),v(train_inds)); %computes Pearson's correlation coefficient
        if frac~=1
            CC_test = corr(double(Xw_test),v(test_inds)); 
        else
            CC_test = nan; 
        end
        
        results.v = v;
        
% % %     CC_train=Xw_train'*v(train_inds)/(norm(Xw_train)*norm(v(train_inds)))
% % %     CC_test=Xw_test'*v(test_inds)/(norm(Xw_test)*norm(v(test_inds)))
% % % end

 
figure; 
plot(y(train_inds),Xw_train,'.b','MarkerSize',25)%,'o'); 
hold all; 
plot(y(test_inds),Xw_test,'.r','MarkerSize',25)%'o'); 
set(gca,'FontSize',18); 
set(gca,'FontName','Times New Roman'); 
ylabel({'Projection score,','in standard deviations'}); 
grid on; 
set(gca,'YTick',[-2*lambda, -lambda, 0, lambda, 2*lambda]);
set(gca,'YTickLabel',[{'-2\sigma','-\sigma','0','\sigma','2\sigma'}]); 

coefficients = polyfit([y(train_inds); y(test_inds)],[Xw_train; Xw_test],1); 

line = coefficients(1).*[y(train_inds); y(test_inds)] + coefficients(2); 
hold on; 
plot([y(train_inds); y(test_inds)], line,'k','LineWidth',5);  %line of best fit computed using all the data

results.wcorr = wcorr;
results.CC_train = CC_train; 
results.CC_test = CC_test; 
results.train_inds = train_inds; 
results.test_inds = test_inds;
y_label = [y(train_inds); y(test_inds)]; projection = [Xw_train; Xw_test];
results.y_label = y_label; results.projection = projection; results.coefficients = coefficients; 

p = 0;
%create testing and training indices
for i = 1:NUM_TESTS
    scrambled = v(randperm(numel(v)));
    train_inds = randsample(num,round(num*frac),'false'); %randomly select 70% of samples for training
    test_inds = 1:num; 
    test_inds(train_inds) = []; 

    %wcorr_perm = Data(:,train_inds)*scrambled(train_inds)/sqrt(scrambled(train_inds)'*Data(:,train_inds)'*Data(:,train_inds)*scrambled(train_inds)); %most correlated direction 
    wcorr_perm = Data(:,train_inds)*(scrambled(train_inds)-mean(scrambled(train_inds)))/sqrt((scrambled(train_inds)-mean(scrambled(train_inds)))'*Data(:,train_inds)'*Data(:,train_inds)*(scrambled(train_inds)-mean(scrambled(train_inds))));

    Xw_train=Data(:,train_inds)'*wcorr_perm; %Projection of the dataset onto the most correlative direction 
    Xw_test=Data(:,test_inds)'*wcorr_perm; 

    %Calculate the correlation coefficient
   
        CC_train_perm = corr(double(Xw_train),scrambled(train_inds)); %computes Pearson's correlation coefficient
        if frac~=1
            CC_test_perm = corr(double(Xw_test),scrambled(test_inds)); 
        else
            CC_test_perm = nan; 
        end
        
% % %     CC_train_perm=Xw_train'*scrambled(train_inds)/(norm(Xw_train)*norm(scrambled(train_inds)));
% % %     CC_test_perm=Xw_test'*scrambled(test_inds)/(norm(Xw_test)*norm(scrambled(test_inds)));
% % %         
    if CC_train_perm > CC_train
        p = p + 1; 
    end
end

p = p/NUM_TESTS;

n = [];
if visualize
    if isequal(vizaxis,'x')
        %n = [-1 -0.5 0 0.5 1]*lambda; %computing the values of independent variable 
        %n = [-1 -0.5 0 0.5 1]*lambda; %computing the values of independent variable
        %n = [mean(y) mean(y)-20].*coefficients(1) + coefficients(2); 
        %since the coefficients are in terms of y, should use y for the x axis scale
        %n = [50.6 58.6 66.5 74.4 82.4].*coefficients(1) + coefficients(2); %if you want to generate images corresponding to independent variable
        %n = [21.1 11.1].*coefficients(1) + coefficients(2); %if you want to generate images corresponding to independent variable
        %n = [mean(v)-2*std(v) mean(v)-std(v) mean(v) mean(v)+std(v) mean(v)+2*std(v)].*coefficients(1) + coefficients(2); %if you want to generate images corresponding to independent variable
        n = [1.9 11.5 21.1 30.7 40.3].*coefficients(1) + coefficients(2);
    else
        n = [-2 -1 0 1 2]*lambda; % if you want to generate images 
    %n = [-3 -2.5 -1.5 0 1.5]*lambda; 
    end
for i = 1:numel(n)
    i
    disp = n(i)*wcorr; %direction in the PCA space. Need to find the direction in the image space
    direction = 0; 
    for j = 1:size(matrix,2)
        load(strcat('eigenvector_',num2str(j))); 
        direction = disp(j)*eigenv + direction; 
    end
    if isequal(mode,'image') || isequal(mode,'Image')
        image{1,i} = imrotate(reshape(direction+A_mean',[M N K]),90);
    else
        %disp = X_mean + n(i)*lambda*wcorr;
        direction = direction + A_mean'; 
        sz = numel(I0);
        
        u = reshape(direction(1:sz),M,N,K);
        v = reshape(direction(sz+1:2*sz),M,N,K); 
        w = reshape(direction(2*sz+1:3*sz),M,N,K); 
        
        f = double(X + u); 
        g = double(Y + v); 
        h = double(Z + w); 

        fields{i,1} = u; 
        fields{i,2} = v; 
        fields{i,3} = w; 

        [dfdx,dfdy,dfdz] = gradient(f); %compute Jacobian map
        [dgdx,dgdy,dgdz] = gradient(g); 
        [dhdx,dhdy,dhdz] = gradient(h);

        %And Jacobian determinant |Du|
        D = (dfdx.*dgdy.*dhdz + dfdy.*dgdz.*dhdx + dfdz.*dgdx.*dhdy - dfdx.*dgdz.*dhdy - dfdy.*dgdx.*dhdz - dgdy.*dhdx.*dfdz); %determinant

        image{1,i} = inpaint_nans3(griddata(double(f),double(g),double(h),double((I0./D)),double(X),double(Y),double(Z))); 
        image{1,i}(image{1,i}>max(I0(:))) = max(I0(:)); 
        image{1,i}(image{1,i}<min(I0(:))) = min(I0(:));
        image{1,i} = imrotate(image{1,i}./sum(image{1,i}(:))*10^6,90); %make sure image is perfectly normalized and rotated
    end
        
end

results.image = image; 
results.fields = fields; %optional output, slows down 

[ B,T] = Gen_Stack(image);
viewdicom(T,0); imcontrast

sz = size(T,2); 
set(gca,'xtick',[sz/10:sz/5:sz])
if isequal(vizaxis,'x')
    n = round((n-coefficients(2))./coefficients(1)*10)/10; %this is to convert back to the y axis proj score 
    %set(gca,'XTickLabel',[{num2str(n(1)), num2str(n(2)), num2str(n(3)), num2str(n(4)), num2str(n(5))}]);
else
    set(gca,'XTickLabel',[{'-2\sigma', '\sigma', '0','\sigma','2\sigma'}]);
    xlabel('Projection score');
end
set(gca,'FontName','Times New Roman');
set(gca,'YTickLabel',[])
ylabel('')
set(gca,'FontSize',18);

%else
    %image = []; n = [];

end

results.n_std = n; 
results.lambda = lambda;
results.p = p; 
results.scores = scores; 
results.covariates = covariates; 
results.matrix = matrix; results.y = y; 

fprintf('r = %.4f, p = %.4f \n', results.CC_train, results.p); 




end

