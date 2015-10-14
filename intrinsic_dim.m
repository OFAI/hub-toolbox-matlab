function no_dims = intrinsic_dim(X)
% Estimate the intrinsic dimensionality of dataset X (Pts x Dims)
%
% Taken from the DR-Toolbox 
% 
% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.7.2b.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, 2010
% University California, San Diego / Delft University of Technology

    X = double(unique(X, 'rows'));
    X = X - repmat(mean(X, 1), [size(X, 1) 1]);
    X = X ./ repmat(var(X, 1) + 1e-7, [size(X, 1) 1]);
    
    % Set neighborhood range to search in
    k1 = 6;
    k2 = 12;

    % Compute matrix of log nearest neighbor distances
    X = X';
    [tmp, n] = size(X);
    X2 = sum(X.^2, 1); 
    knnmatrix = zeros(k2, n);
    if n < 3000
        distance = repmat(X2, n, 1) + repmat(X2', 1, n) - 2 * X' * X;
        distance = sort(distance);
        knnmatrix= .5 * log(distance(2:k2 + 1,:));
    else
        for i=1:n
            distance = sort(repmat(X2(i), 1, n) + X2 - 2 * X(:,i)' * X);
            distance = sort(distance);
            knnmatrix(:,i) = .5 * log(distance(2:k2 + 1))'; 
        end
    end  

    % Compute the ML estimate
    S = cumsum(knnmatrix, 1);
    indexk = repmat((k1:k2)', 1, n);
    dhat = -(indexk - 2) ./ (S(k1:k2,:) - knnmatrix(k1:k2,:) .* indexk);

    % Plot histogram of estimates for all datapoints
    %hist(mean(dhat), 80), pause
    
    % Average over estimates and over values of k
    no_dims = mean(mean(dhat));
end

