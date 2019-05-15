function dataz = pear_fisherz(data)
% Description: transforms pearson R values to Fisher's Z values
%
% input
% -----
% data = vector of pearson R values, 1xN
%
%
% output
% ------
% dataz = vector of transformed fisher's Z values, 1xN 

%% code

sn = sign(data);
dataz = (.5.*log((1+abs(data))./(1-abs(data))));
dataz = dataz .* sn;
