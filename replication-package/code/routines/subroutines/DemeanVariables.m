function [DemeanedXX] = DemeanVariables(X) 
% 
% X is a TxN matrix 

[T,N] = size(X);
MeanX = repmat(mean(X,'omitnan'),T,1);
DemeanedXX = X - MeanX;