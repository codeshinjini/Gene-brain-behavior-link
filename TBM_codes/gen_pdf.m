% Shinjini Kundu (c) 2014
% Transport-Based Morphometry project (TBM)

function [I_out]=gen_pdf(I_in,dc,phi)

%Preprocessing of 3D image

%Inputs:     I_in       3D input image
%            dc         a positive constant
%            phi        kernel for filtering image
%
%Outputs:    I_out      3D output image

I = mat2gray(I_in); %normalize I into a 3D pdf  
I = convn(I,phi,'same'); %faster than imfilter

%I = imfilter(I,phi,'symmetric');%filter the image using kernel phi

%%I = I + dc; 
I = I/sum(I(:)) + dc/numel(I_in); %add constant dc level

I_out = I/sum(I(:)); %normalize I again because it is a 3D pdf    

end


