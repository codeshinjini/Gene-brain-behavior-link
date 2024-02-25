%% Test code to make sure that all folds have projections scores that correspond in direction
for fold = 1:10
    tl = w.train_labels{fold}; 
    pr = w.projected{fold};

    con = (tl ==2); 
    del = (tl ==1); 
    dup = (tl ==3); 
    
    hold all;
    if sum(sign(pr(del,1))) >0  && sum(sign(pr(dup,1)))<0
        fprintf('switch 1\n')
        %dir1 = -dir1;
        pr(:,1) = -pr(:,1); 
    end
    if sum(sign(pr(con,2))) <0 
        fprintf('switch 2 \n'); 
        %dir2 = -dir2;
        pr(:,2) = -pr(:,2);
    end
    
    scatter(pr(del,1),pr(del,2),'yellow'); 
    scatter(pr(con,1),pr(con,2),'magenta')
    scatter(pr(dup,1),pr(dup,2),'cyan')
    legend('del','con','dup')
end

%% Generate visualizations
% Shinjini Kundu (c) 2018, updated 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matter = 'WM'; %or WM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath 'Analysis codes';
clc
k = 3; 

A_mean = w.MEAN;
if isequal(matter,'GM')
    load TBM_gray_down/T1_1_results.mat;
    I0 = results.I0; 
elseif isequal(matter,'WM')
     load TBM_white_down/T1_1_results.mat;
     I0 = results.I0;
end

[M,N,K] = size(I0); 
PLDAprojected = []; 
labels = []; 

for fold = 1:10
    PLDAprojected = [PLDAprojected; w.projected{fold}];
    labels = [labels; w.train_labels{fold}'];
end

[X,Y,Z] = meshgrid(1:N,1:M,1:K); 
sigma1 = std(PLDAprojected(:,1)); 
sigma2 = std(PLDAprojected(:,2)); 

for j = 1:k
    pi(j) = sum(labels==j)/numel(labels);  
    mu_hat(j,:) = sum(PLDAprojected((labels==j),:))./sum(labels==j);
end

V1 = w.dir1; V2 = w.dir2; 

counter = 1;
n1 = [-3 -1.5 0 1.5 3]; n2 = [-3 -1.5 0 1.5 3];
results.n1 = n1; results.n2 = n2; results.sigma1 = sigma1; results.sigma2 = sigma2; 
for i = 1:numel(n1)
    for j = 1:numel(n2)
        fprintf('Now producing the %d th brain \n',counter); 
        direction = A_mean' + n1(i)*sigma1*V1 + n2(j)*sigma2*V2;
        
        index1 = numel(n2)+1-j; index2 = i; %making sure that the images are in the same spatial order as their sampling on the discriminant subspace
        
        for m = 1:k
            distances(m) = norm([n1(i)*sigma1 n2(j)*sigma2]-mu_hat(m,:),2)^2-log(pi(m));
        end
        [~,p] = min(distances); 
        membership(index1,index2) = p; 
             
        %%membership(index1,index2) = predict(mdl,[n1(i)*sigma1 n2(j)*sigma2]); %group membership of each synthetic brain
        
        %visualize the synthetic brain here
        sz = size(direction,1)/3;

        u = reshape(direction(1:sz),M,N,K);
        v = reshape(direction(sz+1:2*sz),M,N,K); 
        z = reshape(direction(2*sz+1:3*sz),M,N,K); 

        f = double(X + u); 
        g = double(Y + v); 
        h = double(Z + z); 

        fields{index1,index2,1} = u; 
        fields{index1,index2,2} = v; 
        fields{index1,index2,3} = z; 

        [dfdx,dfdy,dfdz] = gradient(f); %compute Jacobian map
        [dgdx,dgdy,dgdz] = gradient(g); 
        [dhdx,dhdy,dhdz] = gradient(h);

        %And Jacobian determinant |Du|
        D = (dfdx.*dgdy.*dhdz + dfdy.*dgdz.*dhdx + dfdz.*dgdx.*dhdy - dfdx.*dgdz.*dhdy - dfdy.*dgdx.*dhdz - dgdy.*dhdx.*dfdz); %determinant

        image{index1,index2} = inpaint_nans3(griddata(double(f),double(g),double(h),double((I0./D)),double(X),double(Y),double(Z))); 
        image{index1,index2}(image{index1,index2}>max(I0(:))) = max(I0(:)); 
        image{index1,index2}(image{index1,index2}<min(I0(:))) = min(I0(:)); 
        image{index1,index2} = imrotate(image{index1,index2}./sum(image{index1,index2}(:))*10^6,90);
        counter = counter + 1; 
    end
end

results.image = image; 
results.fields = fields; 
results.membership = membership; 

[B,T] = Gen_Stack(image); 
results.T = T; %images are in the same orer as 
viewdicom(T,0); imcontrast; %the placement of the membership matrix is in the same placement as the images
%y axis corresponds to direction 2

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


%% Analyze displacement fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matter = 'GM'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd Visualization; %load your visualizations
load(strcat(matter,'_TBM_direction.mat')); 

%calculate the scores for standard deviations
sigma1 = std(PROJECTED(:,1)); 
sigma2 = std(PROJECTED(:,2)); 
n = 1.5; %interested in examining the shifts 1.5 standard deviations from the mean

%Calculate the greatest amounts of shift
if isequal(matter,'GM')
    load TBM_gray_down/T1_1_results.mat; 
elseif isequal(matter,'WM')
    load TBM_white_down/T1_1_results.mat; 
end
[M,N,K] =size(results.I0);


mask = results.I0>2*min(results.I0(:)); 
%visualize the synthetic brain here

 for i = 1:2
     if i==1
         direction = direction1; 
     elseif i==2
         direction = direction2; 
     end
     sz = size(direction,1)/3;
     u = reshape(direction(1:sz),M,N,K);
     v = reshape(direction(sz+1:2*sz),M,N,K); 
     z = reshape(direction(2*sz+1:3*sz),M,N,K); 
    
     displacement = (1.5)*(u.^2 + v.^2 + z.^2).^0.5; %1.5 for voxel downsampling

     if i==1
         if isequal(matter,'GM')
            displacement1 = imrotate(displacement(26:end-25,36:end-35,26:end-25),90); 
         elseif isequal(matter,'WM')
            displacement1 = imrotate(displacement(16:end-15,26:end-25,16:end-15),90); 
         end
         displacement1 = displacement1*n*sigma1; %absolute shifts with n std deviation away from mean
         clear displacement;
     elseif i==2
         if isequal(matter,'GM')
            displacement2 = imrotate(displacement(26:end-25,36:end-35,26:end-25),90); 
         elseif isequal(matter,'WM')
            displacement2 = imrotate(displacement(16:end-15,26:end-25,16:end-15),90); 
         end
         displacement2 = displacement2*n*sigma2; %abslute shifts with n std deviation away from mean
         clear displacement
     end
 end
 



%% Calculate visualizations
k = 3;
counter = 1;

if isequal(matter,'GM')
    load TBM_gray_down/T1_1_results.mat;
    I0 = results.I0; 
elseif isequal(matter,'WM')
     load TBM_white_down/T1_1_results.mat;
     I0 = results.I0;
end

[M,N,K] =size(results.I0);


[X,Y,Z] = meshgrid(1:N,1:M,1:K); 
sigma1 = std(PROJECTED(:,1)); 
sigma2 = std(PROJECTED(:,2)); 
n1 = [-3 -1.5 0 1.5 3]; n2 = [-3 -1.5 0 1.5 3];

for j = 1:k
    pi(j) = sum(LABELS==j)/numel(LABELS);  
    mu_hat(j,:) = sum(PROJECTED((LABELS==j),:))./sum(LABELS==j);
end

 for i = 1:numel(n1)
    for j = 1:numel(n2)
        fprintf('Now producing the %d th brain \n',counter); 
        direction = A_mean' + n1(i)*sigma1*direction1 + n2(j)*sigma2*direction2;
        
        index1 = numel(n2)+1-j; index2 = i; %making sure that the images are in the same spatial order as their sampling on the discriminant subspace
        
        for m = 1:k
            distances(m) = norm([n1(i)*sigma1 n2(j)*sigma2]-mu_hat(m,:),2)^2-log(pi(m));
        end
        [~,p] = min(distances); 
        membership(index1,index2) = p; 
             
        %%membership(index1,index2) = predict(mdl,[n1(i)*sigma1 n2(j)*sigma2]); %group membership of each synthetic brain
        
        %visualize the synthetic brain here
        
        sz = size(direction,1)/3;

        u = reshape(direction(1:sz),M,N,K);
        v = reshape(direction(sz+1:2*sz),M,N,K); 
        z = reshape(direction(2*sz+1:3*sz),M,N,K); 

        f = double(X + u); 
        g = double(Y + v); 
        h = double(Z + z); 

        fields{index1,index2,1} = u; 
        fields{index1,index2,2} = v; 
        fields{index1,index2,3} = z; 

        [dfdx,dfdy,dfdz] = gradient(f); %compute Jacobian map
        [dgdx,dgdy,dgdz] = gradient(g); 
        [dhdx,dhdy,dhdz] = gradient(h);

        %And Jacobian determinant |Du|
        D = (dfdx.*dgdy.*dhdz + dfdy.*dgdz.*dhdx + dfdz.*dgdx.*dhdy - dfdx.*dgdz.*dhdy - dfdy.*dgdx.*dhdz - dgdy.*dhdx.*dfdz); %determinant

        image{index1,index2} = inpaint_nans3(griddata(double(f),double(g),double(h),double((I0./D)),double(X),double(Y),double(Z))); 
        image{index1,index2}(image{index1,index2}>max(I0(:))) = max(I0(:)); 
        image{index1,index2}(image{index1,index2}<min(I0(:))) = min(I0(:)); 
        image{index1,index2} = imrotate(image{index1,index2}./sum(image{index1,index2}(:))*10^6,90);
        counter = counter + 1; 
    end
end

 
[B,T] = Gen_Stack(image); 
results.T = T; %images are in the same orer as 
viewdicom(T,0); imcontrast; %the placement of the membership matrix is in the same placement as the images
%y axis corresponds to direction 2
