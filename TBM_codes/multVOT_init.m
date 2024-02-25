%Shinjini Kundu (c) 2022
%wrapper code to run VOT method, updated

function final_results = multVOT_init(I0,I1,f,g,h)
% This function takes an initalization to run
% gradient descent, rather than assuming the identity as the start point.
%% inputs: 
% I0              template image
% I1              source image
% f g h          inital fields to initialize gradient descent. Can leave as
%                  [] if identity is to be used as initialization 
%% outputs: 
% final_results   results of gradient descent with final deformation fields


addpath('DGradient'); 
addpath('codegen/mex/GPExpand');
addpath('codegen/mex/GPReduce'); 

%cd DGradient; 
%mex -O CFLAGS="\$CFLAGS -std=c99" DGradient.c 
%cd ..


numScales = 4; %can set the number of multi-resolution scales for initialization

%1. Initialize potential field of curl-free map (the identity or nearly solution)

[M,N,K]=size(I1);
[X,Y,Z]=meshgrid(1:N,1:M,1:K);

[sx,sy,sz] = size(I0); 

globalIter = 1; 

tic
for scale = numScales:-1:0 
    I0_down = I0; 
    I1_down = I1; 

    f_down = f; 
    g_down = g; 
    h_down = h; 
    
    if scale~=0
        for i = 1:scale
            [X_down1,Y_down1,Z_down1] = meshgrid(1:2^(i-1):sy,1:2^(i-1):sx,1:2^(i-1):sz); 
            I0_down = GPReduce(I0_down); 
            I1_down = GPReduce(I1_down);
            newdim = size(X_down1); %size of the next higher level of the pyramid after I0_down

            if ~isempty(f)
            %downsampled displacement fields
            [X,Y,Z] = meshgrid(1:size(f_down,2),1:size(f_down,1),1:size(f_down,3));
       
            u0_down = 0.5*GPReduce(f_down-X); 
            v0_down = 0.5*GPReduce(g_down-Y); 
            w0_down = 0.5*GPReduce(h_down-Z); 

            %new deformation fields after downsampling
           [X_new,Y_new,Z_new] = meshgrid(1:size(u0_down,2),1:size(u0_down,1),1:size(u0_down,3));
           f_down = u0_down + X_new; 
           g_down = v0_down + Y_new; 
           h_down = w0_down + Z_new; 
            end
        end 
    end
   [X,Y,Z] = meshgrid(1:size(I0_down,2),1:size(I0_down,1),1:size(I0_down,3));


  %initialize the deformation field
    if globalIter==1
        if isempty(f)
            f0 = X; 
            g0 = Y; 
            h0 = Z; 
        else
            %if scale~=0
            %    f0 = u0_down + X; 
            %    g0 = v0_down + Y; 
             %   h0 = w0_down + Z;
            %else 
                f0 = f_down; 
                g0 = g_down; 
                h0 = h_down; 
            %end
        end
   end
   %stop initialization

    if scale > 1
        results = VOT3D_new(I0_down,I1_down,f0,g0,h0,10^-7,scale);
    elseif scale==1
        results = VOT3D_new(I0_down,I1_down,f0,g0,h0,10^-3,scale); 
    elseif scale==0
        results = VOT3D_new(I0_down,I1_down,f0,g0,h0,10^-2,scale); 
    end
    if isempty(results)
        [X,Y,Z] = meshgrid(1:size(I0_down,2),1:size(I0_down,1),1:size(I0_down,3)); 
   %     results.f1 = X; results.f2 = Y; results.f3 = Z; results.I0_recon = I0_down; results.MSE = 10^-3; results.curl = 0; results.I0 = I0_down; results.I1 = I1_down; results.time = 0; 
         results.f1 = f0; results.f2 = g0; results.f3 = h0; 
         [Cx,Cy,Cz] = curl(results.f1,results.f2,results.f3); %compute curl
         [f1x,f1y,f1z]=gradient(results.f1); [f2x,f2y,f2z]=gradient(results.f2);[f3x,f3y,f3z]=gradient(results.f3);detf = (f1x.*f2y.*f3z + f1y.*f2z.*f3x + f1z.*f2x.*f3y - f1x.*f2z.*f3y - f1y.*f2x.*f3z - f1z.*f2y.*f3x);
         results.I0_recon = detf.*interp3(I1_down,results.f1,results.f2,results.f3,'linear',min(I1_down(:)));
  
         %results.I0_recon = I0_down; 
         results.MSE = 10^-3; 
         results.curl = 0.5*(norm(Cx(:),2).^2 + norm(Cy(:),2).^2 + norm(Cz(:),2).^2); results.I0 = I0_down; results.I1 = I1_down; results.time = 0; 
    end
    if globalIter == numScales + 1
        t = toc
        final_results = results;
        final_results.time = t;
        fprintf('the final curl is %d \n', final_results.curl(end)); 
        fprintf('the final MSE is %d \n', final_results.MSE(end)); 
        fprintf('the final time it took was %d minutes \n', t/60); 
        %figure; imagesc(results.I0_recon(:,:,round(K/2))); colorbar; title('morphed'); 
        %figure; imagesc(results.I0(:,:,round(K/2))); colorbar; title('target');
        return; 
    else
        [ ~,~,~,~,~,flag ] = compVOTGradients( results.f1,results.f2,results.f3,zeros(size(I0_down)),zeros(size(I0_down)),0 );
        if (flag)
        %    fprintf('I am not diffeomorphic! \n');
        else
        %    fprintf('I am diffeomorphic \n'); 
        end
        [X2,Y2,Z2] = meshgrid(1:size(X_down1,2),1:size(X_down1,1),1:size(X_down1,3));     
        f0 = 2*GPExpand(results.f1-X,newdim)+X2; 
        g0 = 2*GPExpand(results.f2-Y,newdim)+Y2; 
        h0 = 2*GPExpand(results.f3-Z,newdim)+Z2; 
        [ ~,~,~,~,~,flag ] = compVOTGradients( f0,g0,h0,zeros(size(f0)),zeros(size(f0)),0 );
        if (flag)
            %fprintf('I am not diffeomorphic! \n');
            sigma = 2/pi; 
            [Xt,Yt,Zt]=meshgrid(-3*sigma:3*sigma,-3*sigma:3*sigma,-3*sigma:3*sigma);
            phi = gaussian_bf(Xt,Yt,Zt,sigma); %normalized gaussian kernel in 3D
            f0 = 2*GPExpand(convn(results.f1-X,phi,'same'),newdim)+X2; 
            g0 = 2*GPExpand(convn(results.f2-Y,phi,'same'),newdim)+Y2; 
            h0 = 2*GPExpand(convn(results.f3-Z,phi,'same'),newdim)+Z2;
        else
            %fprintf('I am diffeomorphic \n'); 
        end
        globalIter = globalIter + 1; 
    end
end

end
