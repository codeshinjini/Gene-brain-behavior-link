%Shinjini Kundu (c) 2014
%Computes Cohen's Kappa for contingency tables

function [ kappas] = CohensKappa(cms)
%The function computes Cohen's Kappa, a measure of interrater agreement
%(classifier and gold standard). Because accuracy does not take into account 
%the fact that correct classification could be a result of coincidental
%concordance between the classifiers output and label-generation process.
%The formula is kappa = (P(a) - P(e))/(1-P(e))
%
%Inputs:    cms                     a cell array of confusion matrices 
%
%Outputs:   kappas                  list of kappa values

kappas = zeros(numel(cms),1);

for i = 1:numel(cms)
    confusionmatrix = cms{i};
    dim = size(confusionmatrix,1); 
    colsum = sum(confusionmatrix,1); 
    rowsum = sum(confusionmatrix,2); 
    tot = sum(confusionmatrix(:)); 
    Pa = 0; 
    Pe = 0; 
    for j = 1:dim; 
        Pa = Pa + confusionmatrix(j,j); 
        Pe = Pe + (colsum(j)/tot)*(rowsum(j)/tot);
    end
    Pa = Pa/tot; 
    kappas(i) = (Pa - Pe)/(1-Pe); 
end

