function DIV = new_divergence(f,g,h)
%fourth-order accurate computation of divergence
%inputs:   f g h     deformation fields
%outputs:  DIV       divergence

if isempty(h)
    
[dfdx,~]=new_gradient(f,'fourth'); 
[~,dgdy]=new_gradient(g,'fourth'); 

DIV = dfdx + dgdy; 

else
    
[dfdx,~,~]=new_gradient(f,'fourth'); 
[~,dgdy,~]=new_gradient(g,'fourth'); 
[~,~,dhdz]=new_gradient(h,'fourth');

DIV = dfdx + dgdy + dhdz; 

end


end

