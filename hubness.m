function [Sn, Dk, Nk] = hubness(D, k, isSimilarityMatrix)
% Computes the hubness of a distance matrix using its k nearest neighbors.
% Hubness [1] is the skewness of the n-occurrence histogram.
%
% This file is part of the HUB TOOLBOX available at
% http://ofai.at/research/impml/projects/hubology.html
% (c) 2013, Dominik Schnitzer <dominik.schnitzer@ofai.at>
%
% Usage:
%   Sn = hubness(D) - Computes the hubness (Sk) of the n=5 occurrence histogram
%      (standard)
%
%   [Sn, Dk, Nk] hubness(D, k) - Computes the hubness of the n-occurrence
%      histogram where n (k) is given. Nk is the n-occurrence histogram, Dk
%      are the k nearest neighbors.
%
% [1] Hubs in Space: Popular Nearest Neighbors in High-Dimensional Data
% Radovanovic, Nanopoulos, Ivanovic, Journal of Machine Learning Research 2010

    if (nargin < 2)
        k = 5;
    end
    if (nargin < 3)
        isSimilarityMatrix = 0;
    end
    
    if (isSimilarityMatrix == 1)
        d_self = -Inf;
        sortorder = 'descend';
    else
        d_self = +Inf;
        sortorder = 'ascend';
    end
            
    Dk = zeros(k, size(D, 2));
    for i = 1:size(D, 1)
        % extract distance matrix row
        d = D(i, :);
        
        % make non-finite (NaN, Inf) appear on the end of the sorted list
        d(~isfinite(d)) = d_self;
        d(i) = d_self;
        
        % randomize the distance matrix row to avoid the problem case
        % if all numbers to sort are the same, which would yield high
        % hubness, even if there is none.
        rp = randperm(size(D, 2));
        d2 = d(rp);
        [tmp, d2idx] = sort(d2, sortorder);
        Dk(:, i) = rp(d2idx(1:k));
    end

    Nk = zeros(length(D), 1);
    for i = 1:length(D)
        Nk(i) = sum(sum(Dk == i));
    end

    Sn = local_skewness(Nk);
end

function s = local_skewness(x)
    x0 = x -  mean(x);
    s2 = mean(x0.^2);
    m3 = mean(x0.^3);
    s = m3 ./ s2.^(1.5);
end
