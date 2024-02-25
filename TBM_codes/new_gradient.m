%Computes gradient using second-order methods
%Shinjini Kundu (c) 2015

function [dIdx, dIdy, dIdz] =  new_gradient(image, option)
%inputs:  image              a 3D image
%         option             'first' or 'second'
%                            choosing 'first' is the same as the result
%                            obtained by matlab's gradient function
%output:  dIdx, dIdy, dIdz   gradients in x,y and z directions

%note: assumes that image has same size in all 3 dimensions

dim = numel(size(image));

if iscolumn(image)||isrow(image)
    dim =1; 
end

if dim==1
    if ~isrow(image)
        image = image'; 
    end
    M = numel(image); 
elseif dim==2
    [M,N] = size(image);
    dIdx = zeros(size(image)); 
    dIdy = zeros(size(image));
elseif dim==3
    [M,N,K] = size(image); 
    dIdx = zeros(size(image)); 
    dIdy = zeros(size(image)); 
    dIdz = zeros(size(image)); 
end

if isequal(option, 'first')
    %centered difference, first order gradient approximation
    G = sparse(0.5*full(gallery('tridiag',M,1,0,-1)));
    G(:,1) = [-1; 1; zeros(M-2,1)];
    G(:,M) = [zeros(M-2,1); -1; 1]; 

    if dim==1
        dIdx = (image*G)'; 
    elseif dim==2
        dIdx = image*G; 
        dIdy = (image'*G)'; 
        dIdz = [];
    elseif dim==3
        for i = 1:K
            dIdx(:,:,i) = image(:,:,i)*G; 
            dIdy(:,:,i) = (image(:,:,i)'*G)'; 
            dIdz(:,i,:) = (squeeze(image(:,i,:))*G); 
        end
    end
end

if isequal(option,'second')
    %second order gradient approximation, 

    g = sparse((full(gallery('circul',[1; -8; 0; 8; -1; zeros(M-5,1)]))./12)');
    G(:,1) = [-1; 1; zeros(M-2,1)];
    G(:,2) = [-0.5; 0; 0.5; zeros(M-3,1)]; 
    G(:,3:M-2) = g(:,1:M-4); 
    G(:,M-1) = [zeros(M-3,1); -0.5; 0; 0.5]; 
    G(:,M) = [zeros(M-2,1); -1; 1]; 

    if dim==1
        dIdx = (image*G)'; 
    elseif dim==2
        dIdx = image*G; 
        dIdy = (image'*G)'; 
        dIdz = [];
    elseif dim==3
        for i = 1:K
            dIdx(:,:,i) = image(:,:,i)*G; 
            dIdy(:,:,i) = (image(:,:,i)'*G)'; 
            dIdz(:,i,:) = (squeeze(image(:,i,:))*G); 
        end
    end
    
end

if isequal(option,'third')
    %second order gradient approximation, 

    g = sparse((full(gallery('circul',[-1/60; 3/20; -3/4; 0; 3/4; -3/20; 1/60; zeros(M-7,1)])))');
    G(:,1) = [-1; 1; zeros(M-2,1)];
    G(:,2) = [-0.5; 0; 0.5; zeros(M-3,1)]; 
    G(:,3) = [1/12; -2/3; 0; 2/3; -1/12; zeros(M-5,1)]; 
    G(:,4:M-3) = g(:,1:M-6);
    G(:,M-2) = [zeros(M-5,1); 1/12; -2/3; 0; 2/3; -1/12];  
    G(:,M-1) = [zeros(M-3,1); -0.5; 0; 0.5]; 
    G(:,M) = [zeros(M-2,1); -1; 1]; 


    if dim==1
        dIdx = (image*G)'; 
    elseif dim==2
        dIdx = image*G; 
        dIdy = (image'*G)';
        dIdz = [];
    elseif dim==3
        for i = 1:K
            dIdx(:,:,i) = image(:,:,i)*G; 
            dIdy(:,:,i) = (image(:,:,i)'*G)'; 
            dIdz(:,i,:) = (squeeze(image(:,i,:))*G); 
        end
    end
end

if isequal(option,'fourth')
    %second order gradient approximation, 

    g = sparse((full(gallery('circul',[1/280; -4/105; 1/5; -4/5; 0; 4/5; -1/5; 4/105; -1/280; zeros(M-9,1)])))');
    G(:,1) = [-1; 1; zeros(M-2,1)];
    G(:,2) = [-0.5; 0; 0.5; zeros(M-3,1)]; 
    G(:,3) = [1/12; -2/3; 0; 2/3; -1/12; zeros(M-5,1)];
    G(:,4) = [-1/60; 3/20; -3/4; 0; 3/4; -3/20; 1/60; zeros(M-7,1)]; 
    G(:,5:M-4) = g(:,1:M-8);
    G(:,M-3) = [zeros(M-7,1); -1/60; 3/20; -3/4; 0; 3/4; -3/20; 1/60]; 
    G(:,M-2) = [zeros(M-5,1); 1/12; -2/3; 0; 2/3; -1/12];  
    G(:,M-1) = [zeros(M-3,1); -0.5; 0; 0.5]; 
    G(:,M) = [zeros(M-2,1); -1; 1]; 


    if dim==1
        dIdx = (image*G)'; 
    elseif dim==2
        dIdx = image*G; 
        dIdy = (image'*G)';
        dIdz = [];
    elseif dim==3
        for i = 1:K
            dIdx(:,:,i) = image(:,:,i)*G; 
            dIdy(:,:,i) = (image(:,:,i)'*G)'; 
            dIdz(:,i,:) = (squeeze(image(:,i,:))*G); 
        end
    end
end




end

