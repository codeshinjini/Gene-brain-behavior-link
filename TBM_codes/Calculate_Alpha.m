function [Thresh,error_subspace]=Calculate_Alpha(FINAL_FEATS,labels,nPLDA)
%% This piece of code runs the PLDA function with different 
%  values of alpha and plot two error curves
%  The projection metric distance of two consequent subspaces
%  For    
%% Load the data set
[Npnt,Nfeatures] = size(FINAL_FEATS);
%%

x = 2.^(-10:10); 
counter = 0; 

for i = 1:Npnt
    fprintf('Now loading eigenvector %d \n', i); 
    load(strcat('eigenvector_',num2str(i)));
    VecPCA(:,i) = eigenv; 
end

fprintf('Done load the eigenvectors, now computing alpha values \n'); 
for Alpha=x
    Alpha
    counter=counter+1;
    [PLDA_directions,~] = PLDA(FINAL_FEATS', labels, Alpha, nPLDA); 
    Vec(:,:,counter) = VecPCA*PLDA_directions; 
    if counter>1
          Vec(:,:,counter)=Vec(:,:,counter)*diag(sign(diag(Vec(:,:,counter)'*Vec(:,:,counter-1))));               
          error_subspace(counter-1)= Projection_metric( Vec(:,:,counter),Vec(:,:,counter-1)); 
    end
end 
%% Calculate twice the half life of Alpha,

%Fit an exponential to log(error_subspace)
f =@(a,b,x) log(a)-b*x; 
options = fitoptions('Method','linearLeastSquares');
F_fitted = fit(x(1:end-1)',log(error_subspace)', f, ...
    'StartPoint', [1,x(1)], ...
    'Lower', [0,0],'Robust','LAR');
coeff=coeffvalues(F_fitted);%Get coefficients of fitted function
Thresh=log(2)*(1/coeff(2));%Calculate twice the half life = 2(log(2)/b)

figure
plot(x(1:end-1),error_subspace,'linewidth',2)
title('Stability of subspace w.r.t. Alpha','fontsize',24)
ylabel('Projection metric between two consequent subspaces','fontsize',20)
xlabel('Alpha','fontsize',20)
set(gca,'fontsize',20)
yL = get(gca,'YLim');
hold on
line([Thresh Thresh],yL,'Color','r','linewidth',2);
grid on
