function [A,B,C,determinant] = IOT3Dnew(I0,I1)
%This function receives I0 and I1 as 3D input images and find an initial
%masspreserving mapping between them. The details for the method can be
%found in,
%
%  Haker, Steven, et al. "Optimal mass transport for registration and warping." 
%  International Journal of Computer Vision 60.3 (2004): 225-240.
%
%The mass preserving mapping is from I0 to I1 (Hence it pushes I1 to I0) and 
%provides the following equation
%
%   |D_{A,B,C}(x,y,z)|I1(A(x),B(x,y),C(x,y,z))=I0(x,y,z),
%
%where D_{A,B,C} is the Jacobian matrix of the mapping
%f=[A(x),B(x,y),C(x,y,z)] and |.| denotes the determinant operator. 
%
%Inputs:
%             I0 :      3D source image
%             I1 :      3D target image
%Outputs:     
%             A               transport plan in x dimension
%             B               transport plan in y dimension
%             C               transport plan in z dimension
%             determinant     determinant of the Jacobian matrix
%
%Written by Soheil Kolouri � 2014

close all

[M,N,K] = size(I0);
x = 1:N;  %1-Dimensional x, Domain of J0 and J1.
y = 1:M;  %1-Dimensional y
z = 1:K;  %1-Dimensional z

[X,Y,Z] = meshgrid(x,y,z);%Grid of the I0 and I1.

%Integrate over z domain, resulting in a single 2D images in the x-y plane 

if size(I0)~=size(I1)
    error('dimension mismatch'); 
end

% Sum over z to get 2D images i0 and i1
i0 = sum(I0,3); 
i1 = sum(I1,3);

% Sum over y to get 1D signals J0 and J1

J0 = sum(i0,1);
J1 = sum(i1,1);

clear i1

fprintf('Now initial 1D signals are found! \n'); 

%%Then we find the mass preserving mapping from J0 to J1, 
%                 a'(x)J1(a(x))=J0(x)
% where a' is the gradient of a.

[a,da] = IOT1D(J0,J1);% compute a(x) function

fprintf('Computed 1D transform for the initial signals! \n'); 

% We form A(x,y,z)=a(x), and dA to be its 

A = repmat(single(a),[M,1,K]);       
dA =repmat(single(da),[M,1,K]);      

clear a da

IA = dA.*interp3(single(I1),single(A),single(Y),single(Z),'spline');%Interpolate I1 and call it IA
iA = sum(IA,3); %Take the sum of IA with respect to Z and name it iA

clear dA

fprintf('Now done projecting the set of 1D signals back \n'); 

%%Then we find b(a(x),y) for all y.
b = zeros(M,N); 

% Compute b(a(x),y) portions

fprintf('Starting to compute the set of 1D signals for the entire 3D image \n'); 

%figure
parfor i=1:N    
    J0 = i0(:,i);%Let J0 be the i'th column of I0 
    J1 = iA(:,i);%Let J1 be the i'th column of IA (Look at line 68)
    [b(:,i),db(:,i)] = IOT1D(J0,J1);                             
end

clear i0 J1 J0 iA

fprintf('Done computing the set of 1D signals for the entire 3D image \n'); 

% Generate B(A(x,y,z),y,z)=b(a(x),y) and its jacobian matrix

B = repmat(single(b),[1,1,K]);       %A is [a(1) ... a(N);...;a(1)... a(N)] 
dB =repmat(single(db),[1,1,K]);      %A is [da(1) ... da(N);...;da(1)... da(N)] 

clear b db

IAB = dB.*ba_interp3(double(IA),double(X),double(B),double(Z),'cubic');%Interpolate IA and call it IAB

clear IA dB

fprintf('Done interpolating back to the 3D image \n'); 

% Show results so far
% IABt=IAB/max(IAB(:));
% for i=1:64;
%      imshow(IABt(:,:,i));
%      pause(.1)
%  end

% Compute c(a(x),b(x,y),z) portions 
for i = 1:M
    for j = 1:N 
         [C(i,j,:),dC(i,j,:)] = IOT1D(squeeze(I0(i,j,:)),squeeze(IAB(i,j,:)));                  
    end
end

clear K M N IAB dC

%%% Uncomment to check the results
% IABC = dC.*interp3(IAB,X,Y,C,'spline');%Interpolate I1 and call it IA
% IABCt=IABC/max(IABC(:));
% for i=1:64;
%      imshow(IABCt(:,:,i));
%      pause(.1)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Up to here we have calculated
% 
%       [A(X,Y,Z),B(A(X,Y,Z),Y,Z), C(A(X,Y,Z),B(A(X,Y,Z),Y,Z),Z)]
%
% Hence we need to regrid the mappings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Now regridding the images'); 

C = griddata(double(A),double(B),double(Z),double(C),double(X),double(Y),double(Z),'linear'); %Regrid c to get C

B = griddata(double(A),double(Y),double(Z),double(B),double(X),double(Y),double(Z),'linear'); %Regrid c to get C

clear X Y Z

fprintf('Done regridding. Now inpainting \n'); 

if sum(isnan(C(:)))
    C=inpaint_nans3(C);
end

if sum(isnan(B(:)))
    B=inpaint_nans3(B);
end

fprintf('Now compute Jacobian determinant \n'); 
% Now we calculate the determinant of the Jacobian matrix for the regridded 
% mapping.

[dAdx,dAdy,dAdz] = gradient(A);%get gradient of A 
[dBdx,dBdy,dBdz] = gradient(B);%get gradient of B
[dCdx,dCdy,dCdz] = gradient(C);%get gradient of C

determinant = dAdx.*dBdy.*dCdz + dAdy.*dBdz.*dCdx + dAdz.*dBdx.*dCdy - dAdz.*dBdy.*dCdx - dAdy.*dBdx.*dCdz - dAdx.*dBdz.*dCdy;
clear dAdx dAdy dAdz dBb dBdz dBdy dBdz dCdx dCdy dCdz

%%%% Uncomment to check the final result.
Iout = determinant.*interp3(single(I1),single(A),single(B),single(C),'spline'); %Interpolate I1
% Ioutt= Iout/max(Iout(:));
% figure
% for i=1:K
%     imshow(Ioutt(:,:,i))
%     pause(.3)
% end 

Ierror= mean((Iout(:)/sum(Iout(:))-I0(:)/sum(I0(:))).^2);
['The MSE between |Df|I1(f) and I0 is equal to, error= ',num2str(Ierror)]
clear Iout Ierror I0 I1

end
