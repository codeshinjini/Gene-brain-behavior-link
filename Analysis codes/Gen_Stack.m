%Shinjini Kundu (c) 2015
function [ B,T ] = Gen_Stack( IMAGES )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

slices = size(IMAGES{1,1},3);

for k = 1:slices
    BIG = []; 
    MEAN = []; 
    for i =  1:size(IMAGES,1)
        BIG1 = []; MEAN1 = []; 
        center_image = IMAGES{i,round(size(IMAGES,2)/2)}; 
        for j = 1:size(IMAGES,2)
             BIG1 = [BIG1 IMAGES{i,j}(20:end-20,20:end-20,k)-center_image(20:end-20,20:end-20,k)];
             MEAN1 = [MEAN1 center_image(20:end-20,20:end-20,k)];

            %BIG1 = [BIG1 IMAGES{i,j}(3:end-3,140:end-140,k)-center_image(3:end-3,140:end-140,k)];
            %MEAN1 = [MEAN1 center_image(3:end-3,140:end-140,k)];
            % BIG1 = [BIG1 IMAGES{i,j}(3:end-3,120:end-120,k)-center_image(3:end-3,120:end-120,k)];
            % MEAN1 = [MEAN1 center_image(3:end-3,120:end-120,k)]; 
        end
        BIG = [BIG; BIG1];
        MEAN = [MEAN; MEAN1]; 
    end

    B(:,:,k) = BIG; 
    %M(:,:,k) = MEAN; 
    T(:,:,k) = BIG + MEAN; 
end

T = mat2gray(T); 
B = mat2gray(B); 

end

