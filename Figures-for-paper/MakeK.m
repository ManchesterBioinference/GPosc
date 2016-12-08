function [ K ] = MakeK( cov,x )
%MAKEK makes covariance matrix

CovMatrix1 = zeros(length(x),length(x));

for i = 1:length(x)
    CovMatrix1(:,i) = circshift(cov',i-1);

end
    CovMatrix1 = CovMatrix1';
    
K = triu(CovMatrix1,0)+(triu(CovMatrix1,1))';


end

