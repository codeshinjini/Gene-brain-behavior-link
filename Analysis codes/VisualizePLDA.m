%Visualization code 
% Shinjini Kundu (c) 2018

function results = VisualizePLDA(FINAL_FEATS, labels,k,I0,A_mean)

nPLDA= 2; 
alpha = 1;
type = 'brain'; 

if min(labels)~=1
    labels = 1-min(labels)+ labels;
end

[M,N,K] = size(I0); 
[Vec, ~] = PLDA(FINAL_FEATS', labels, alpha, nPLDA);
PLDAprojected = FINAL_FEATS*Vec;

%compute the sample synthetic brains
fprintf('Now computing sample synthetic brains \n'); 
v1 = Vec(:,1);%x direction
v2 = Vec(:,2);%y direction

V1 = zeros(3*M*N*K,1); V2 = V1; 
for i = 1:size(FINAL_FEATS,2)
    load(strcat('eigenvector_',num2str(i))); 
    V1 = V1 + v1(i)*eigenv; 
    V2 = V2 + v2(i)*eigenv;
end

counter = 1;
%mdl = fitcknn(PLDAprojected,labels,'NumNeighbors',k);
[X,Y,Z] = meshgrid(1:N,1:M,1:K); 
sigma1 = std(PLDAprojected(:,1)); 
sigma2 = std(PLDAprojected(:,2)); 

for j = 1:k
    pi(j) = sum(labels==j)/numel(labels);  
    mu(j,:) = sum(FINAL_FEATS((labels==j),:))./sum(labels==j);
    mu_hat(j,:) = mu(j,:)*Vec;
end

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
        w = reshape(direction(2*sz+1:3*sz),M,N,K); 

        f = double(X + u); 
        g = double(Y + v); 
        h = double(Z + w); 

        fields{index1,index2,1} = u; 
        fields{index1,index2,2} = v; 
        fields{index1,index2,3} = w; 

        [dfdx,dfdy,dfdz] = gradient(f); %compute Jacobian map
        [dgdx,dgdy,dgdz] = gradient(g); 
        [dhdx,dhdy,dhdz] = gradient(h);

        %And Jacobian determinant |Du|
        D = (dfdx.*dgdy.*dhdz + dfdy.*dgdz.*dhdx + dfdz.*dgdx.*dhdy - dfdx.*dgdz.*dhdy - dfdy.*dgdx.*dhdz - dgdy.*dhdx.*dfdz); %determinant

        image{index1,index2} = inpaint_nans3(griddata(double(f),double(g),double(h),double((I0./D)),double(X),double(Y),double(Z))); 
        image{index1,index2}(image{index1,index2}>max(I0(:))) = max(I0(:)); 
        image{index1,index2}(image{index1,index2}<min(I0(:))) = min(I0(:)); 
        if isequal(type,'brain')
            image{index1,index2} = imrotate(image{index1,index2}./sum(image{index1,index2}(:))*10^6,90);
        else
            image{index1,index2} = image{index1,index2}./sum(image{index1,index2}(:))*10^6; %make sure image is perfectly normalized 
        end
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

end
