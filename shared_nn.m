function Dsnn = shared_nn(D, k)
% Transforms the given distance matrix into new one using a shared nearest
% neighbor transform with the given neighborhood radius k.
%
% This file is part of the HUB TOOLBOX available at
% http://ofai.at/research/impml/projects/hubology.html
% (c) 2013, Dominik Schnitzer <dominik.schnitzer@ofai.at>
%
% Usage:
%   Dsnn = shared_nn(D, k) - Use SNN with a neighborhood radius of k on the
%      given distance matrix. The new distances are returned in Dsnn.

    if (nargin < 2)
        k = 10;
        fprintf('No neighborhood radius given. Using k=10\n');
    end
    
    n = size(D, 1);
    z = false(size(D));
    for i = 1:n
        di = D(i, :);
        di(i) = +Inf;
        [tmp, nn] = sort(di);
        z(i, nn(1:k)) = 1;
    end
    
    Dsnn = zeros(size(D), class(D));
    for i = 1:(n-1)
        zi = z(i,:);
        j_idx = i+1:n;
        
        % BSXFUN
        % Dij = sum(bsxfun(@and, zi, z(j_idx, :)), 2);
        
        % REPMAT
        Dij = sum(repmat(zi, length(j_idx), 1) & z(j_idx, :), 2);
        
        Dsnn(i, j_idx) = 1 - Dij./k;
        Dsnn(j_idx, i) = Dsnn(i, j_idx);
    end
    
end
