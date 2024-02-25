function [ dx,dy,dz ] = cgradient( f )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[dx, dy, dz] = gradient(f); 

dx(:,1,:) = floor(dx(:,1,:)); 
dx(:,end,:) = ceil(dx(:,end,:)); 

dy(1,:,:) = floor(dy(1,:,:)); 
dy(end,:,:) = ceil(dy(end,:,:)); 

dz(:,:,1) = floor(dz(:,:,1)); 
dz(:,:,end) = ceil(dz(:,:,end)); 

end

