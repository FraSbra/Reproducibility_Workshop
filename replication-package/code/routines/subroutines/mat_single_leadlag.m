function [LeadLagX] = mat_single_leadlag(X,leadlagOrder)

[T,N]=size(X); 

LeadLagX = NaN(T,N);

if leadlagOrder>=0
    LeadLagX(1:end-leadlagOrder,:) = X(leadlagOrder+1:end,:);
elseif leadlagOrder<0
    lagOrder = abs(leadlagOrder);
    LeadLagX(lagOrder+1:end,:) = X(1:end-lagOrder,:);
end
