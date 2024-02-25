%Shinjini Kundu (c) 2022
%Kundu et al approach to solving for mass-preserving mapping, updated

function [results] = VOT3D_new(I0,I1,f1,f2,f3,step_size,scale)

addpath('freezeColors')
%addpath F:/TBM_codes/freezeColors;

iter = 1;
lambda =3; %10; %51200; %50; %3.5*10^-12; 

cutoff = 10^-3; %2*10^-4; %4*10^-4;
if scale==0
    cutoff = 3.5*10^-3; %5*10^-3
end

if isempty(step_size)
    step_size = 10^-4; %10^12; 
end
converged = 0;

results = []; 
sigma = 1;
DC_level = 0.1; %0.1; %Create Gaussian kernel in 3D

% % %use this code for 2015 matlab version and beyond
% I0 = mat2gray(imgaussfilt3(I0,sigma)) + DC_level; 
% I1 = mat2gray(imgaussfilt3(I1,sigma)) + DC_level;

%else, can use this 2014 matlab version
[Xt,Yt,Zt]=meshgrid(-3*sigma:3*sigma,-3*sigma:3*sigma,-3*sigma:3*sigma);
phi = gaussian_bf(Xt,Yt,Zt,sigma); %normalized gaussian kernel in 3D
I0 = mat2gray(convn(I0,phi,'same')) + DC_level; 
I1 = mat2gray(convn(I1,phi,'same')) + DC_level; 

I0 = I0./sum(I0(:)); %can normalize to 10^6
I1 = I1./sum(I1(:)); 

[M,N,K]=size(I1);
[X,Y,Z]=meshgrid(1:N,1:M,1:K);
figure(1)

mask = ones(size(I0)); 
for i = 1:size(I0,3)
    mask(:,:,i) = im2bw(I0(:,:,i),min(I0(:)));
end

I0 = I0*10^6; 
I1 = I1*10^6;

while(true)
    %fprintf('Now on interation %d \n', iter); 
    if iter ==1
        [ f1t,f2t,f3t,I0_recon,Ierror,flag ] = compVOTGradients( f1,f2,f3,I0,I1,lambda ); 
        if (flag)
            error('The initial deformation field is not diffeomorphic'); 
        else
            xk1_temp = f1-step_size*f1t; xk2_temp = f2-step_size*f2t; xk3_temp = f3-step_size*f3t;
            yk1_temp = xk1_temp; yk2_temp = xk2_temp; yk3_temp = xk3_temp;
            
            %check to make sure that the updated fields are diffeomorphic
            [~,~,~,~,~,flag] = compVOTGradients(yk1_temp,yk2_temp,yk3_temp,I0,I1,lambda);
            %if not diffeomorphic, need to take a smaller stepsize
            while(flag && ~converged)
                step_size = step_size/2; 
                if step_size < (10^-8) %if there is no stepsize that will enable a diffeomorphic deformation, you have converged
                    converged = 1;
                    step_size = 0;
                    results.f1 = f1; 
                    results.f2 = f2; 
                    results.f3 = f3; 
                end
                xk1_temp = f1-step_size*f1t; xk2_temp = f2-step_size*f2t; xk3_temp = f3-step_size*f3t;
                yk1_temp = xk1_temp; yk2_temp = xk2_temp; yk3_temp = xk3_temp;
                [~,~,~,~,~,flag] = compVOTGradients(yk1_temp,yk2_temp,yk3_temp,I0,I1,lambda);
            end
            xk1 = xk1_temp; xk2 = xk2_temp; xk3 = xk3_temp;
            yk1 = yk1_temp; yk2 = yk2_temp; yk3 = yk3_temp; 
            %fprintf('the stepsize is %d \n', step_size);
               
            yk1minus1 = zeros(size(I0)); yk2minus1 = zeros(size(I0)); yk3minus1 = zeros(size(I0)); 
            xk1minus1 = zeros(size(I0)); xk2minus1 = zeros(size(I0)); xk3minus1 = zeros(size(I0)); 
        end
    end
   
    
    if iter > 1
        [ f1t,f2t,f3t,I0_recon,Ierror ] = compVOTGradients( yk1minus1,yk2minus1,yk3minus1,I0,I1,lambda ); 
        xk1_temp = yk1minus1-step_size*f1t; xk2_temp = yk2minus1-step_size*f2t; xk3_temp = yk3minus1-step_size*f3t;
        yk1_temp = xk1_temp + (iter-2)/(iter+1)*(xk1_temp - xk1minus1); yk2_temp = xk2_temp + (iter-2)/(iter+1)*(xk2_temp - xk2minus1); yk3_temp = xk3_temp + (iter-2)/(iter+1)*(xk3_temp - xk3minus1);
        
        %check whether updated fields are diffeomorphic
        [~,~,~,~,~,flag] = compVOTGradients(yk1_temp,yk2_temp,yk3_temp,I0,I1,lambda);
        while (flag && ~converged)
            step_size = step_size/2; 
            if step_size < (10^-8) %if there is no stepsize that will enable a diffeomorphic deformation, you have converged
                converged = 1; 
                step_size = 0; 
                if iter==2
                    results.f1 = f1; 
                    results.f2 = f2; 
                    results.f3 = f3;
                end
            end
            xk1_temp = yk1minus1-step_size*f1t; xk2_temp = yk2minus1-step_size*f2t; xk3_temp = yk3minus1-step_size*f3t;
            yk1_temp = xk1_temp + (iter-2)/(iter+1)*(xk1_temp - xk1minus1); yk2_temp = xk2_temp + (iter-2)/(iter+1)*(xk2_temp - xk2minus1); yk3_temp = xk3_temp + (iter-2)/(iter+1)*(xk3_temp - xk3minus1);
            [~,~,~,~,~,flag] = compVOTGradients(yk1_temp,yk2_temp,yk3_temp,I0,I1,lambda);
        end
        xk1 = xk1_temp; xk2 = xk2_temp; xk3 = xk3_temp;
        yk1 = yk1_temp; yk2 = yk2_temp; yk3 = yk3_temp; 
       % fprintf('the stepsize is %d \n', step_size); 
    end
    
    if (~converged)
        yk1minus2 = yk1minus1; yk2minus2 = yk2minus1; yk3minus2 = yk3minus1; 
        yk1minus1 = yk1; yk2minus1 = yk2; yk3minus1 = yk3; 
        xk1minus1 = xk1; xk2minus1 = xk2; xk3minus1 = xk3; 
    end
    
    %plotting code
    subplot(221)
    imshow(squeeze(sum(I0_recon,1)),[]); 
    title('$$det(D{\bf f})I_1({\bf f})$$','interpreter','latex','fontsize',20)
    freezeColors
    subplot(222)
    if iter==1
        showgrid(squeeze(X(:,:,round(K/2))-f1(:,:,round(K/2))),squeeze(Y(:,:,round(K/2))-f2(:,:,round(K/2))),3)
    else
        showgrid(squeeze(X(:,:,round(K/2))-yk1minus2(:,:,round(K/2))),squeeze(Y(:,:,round(K/2))-yk2minus2(:,:,round(K/2))),3)
    end
    title('$${\bf x}-{\bf f}({\bf x})$$','interpreter','latex','fontsize',20)
    subplot(2,2,3)
    %err(iter)=.5*sum(((Ierror(:)./I0(:))).^2)/numel(I0(:));
    err(iter)=.5*sum(((Ierror(:)./I0(:)).*mask(:)).^2)/nnz(mask(:)); %numel(I0(:)); %relative MSE
    plot(err,'linewidth',2)
    title('MSE: $$\frac{1}{2}\|det(D{\bf f})I_1({\bf f})-I_0\|^2$$','interpreter','latex','fontsize',20)
    grid on
    subplot(2,2,4)
    if iter==1
        [Cx,Cy,Cz,~] = curl(f1,f2,f3); 
    else
        [Cx,Cy,Cz,~] = curl(yk1minus2,yk2minus2,yk3minus2); 
    end
    errorcurl(iter)=0.5*(norm(Cx(:),2).^2 + norm(Cy(:),2).^2 + norm(Cz(:),2).^2);
    plot(errorcurl,'r','linewidth',2)
    title('Curl: $$\frac{1}{2}\|\nabla\times {\bf f}\|^2$$','interpreter','latex','fontsize',20)
    grid on
    drawnow
% %     %end of plotting code
    
    if (converged || err(iter) <= cutoff)
        return;
    end
    results.f1 = yk1minus2; 
    results.f2 = yk2minus2; 
    results.f3 = yk3minus2; 
    results.I0_recon = I0_recon; 
    results.MSE = err; 
    results.curl = errorcurl; 
    results.I0 = I0; 
    results.I1 = I1;
    results.cutoff = cutoff; 
    results.lambda = lambda;
    results.sigma = sigma; 
    results.DC_level = DC_level;
    results.sigma = sigma;
    iter = iter + 1;
    

    
end
    

end 


%     [ f1t,f2t,f3t,I0_recon,Ierror ] = compVOTGradients( f1,f2,f3,I0,I1 ); 
%     
% %    step_size=.05/(max(sqrt(f1t(:).^2+f2t(:).^2+f3t(:).^2)));
%     f1=f1-step_size*f1t;
%     f2=f2-step_size*f2t;
%     f3=f3-step_size*f3t; 
      
