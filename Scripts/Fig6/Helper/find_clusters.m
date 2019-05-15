function varargout = find_clusters(x)

% Will find clusters of 0 and 1 in a time-series and report as numeric
% output from 1 to n_cluster
%
% Output 1: Cluster sizes
% Output 2: Cluster start index
% Output 3: Cluster end index
% Output 4: Clusters labeled
%
% Example:
% x = [0 1 1 1 0 0 1 1 0 0];
% clust_size = find_clusters(x);
% clust_size =
%
%     3     2
%
% [clust_size,first,last] = find_clusters(x);
% start =
% 
%      2     7
% 
% 
% stop =
% 
%      4     8
%
% [clust_size,start,stop,clust_labels] = find_clusters(x)
% clust_labels =
% 
%      0     1     1     1     0     0     2     2     0     0

% TODO: generalize to matrix (add row of 0 to the end to disconnect
% clusters between different iterations, in for loop reset index between columns using mod)

flipon = 0;
if size(x,2)==1
    flipon = 1;
    x = x';
end

% for now clustering for dummies
varargout{2} = find(diff([0 x 0])==1); % cluster onset
varargout{3} = find(diff([0 x 0])==-1)-1; % cluster offset

% if varargout{3}(1)==1
%     varargout{2} = [1 varargout{2}];
% end

varargout{1} = varargout{3}-varargout{2}+1; % cluster size
if isempty(varargout{1})
    varargout{1} = 0; % no clusters
end
if flipon
    for i = 1:3
        varargout{i} = varargout{i}';
    end
end

if nargout <= 3
    return
end

varargout{4} = zeros(1,length(x));
for i_clust = 1:length(varargout{2})
    varargout{4}(varargout{2}(i_clust):varargout{3}(i_clust)) = i_clust;
    if flipon
        varargout{4} = varargout{4}';
    end
end