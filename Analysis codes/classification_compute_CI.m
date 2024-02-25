%Shinjini Kundu (c) 2022

function [average, lowerCI, upperCI] = classification_compute_CI(vector,T)
%computes confidence interval for a vector assuming normal approximation of the binomal
%inputs:     vector     the vector of values 
%outputs:    mean       mean
%            lowerCI    lower range of the 95% confidence interval
%            upperCI    upper range of the 95% confidence interval

average = mean(vector); 

error = 1-mean(vector); 

range = 1.96*sqrt(error*(1-error)/T); 
lowerCI = round((1-(error+range))*10^3)/10^3; 
upperCI = round((1-(error-range))*10^3)/10^3; 

end

