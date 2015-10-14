function Dmp = mutual_proximity(D, type)
% Applies Mutual Proximity (MP) [1] on a distance matrix. The return value is
% converted to a distance matrix again. The resulting distance matrix
% should show lower hubness.
%
% This file is part of the HUB TOOLBOX available at
% http://ofai.at/research/impml/projects/hubology.html
% (c) 2013, Dominik Schnitzer <dominik.schnitzer@ofai.at>
%
% Usage:
%   Dmp = mutual_proximity(D, type) - Applies MP on the distance matrix 'D'
%      using the selected variant ('type'). The transformed distance matrix
%      is returned.
%
% Possible types:
%   'empiric': Uses the Empirical distribution to perform Mutual Proximity.
%   'gauss': (requires the Statistics Toolbox (the mvncdf() function)
%      Assumes that the distances are Gaussian distributed.
%   'gaussi': Assumes that the distances are independently Gaussian
%      distributed. (fastest Variante)
%   'gammai': Assumes that the distances follow a Gamma distribution and
%      are independently distributed.
%
% [1] Local and global scaling reduce hubs in space, 
% Schnitzer, Flexer, Schedl, Widmer, Journal of Machine Learning Research 2012

    if (nargin < 2)
        mp_func = @mp_empiric;
        fprintf('No Mutual Proximity type given. Using: ''empiric''\n');
        fprintf('For fast results use: ''gaussi''\n');
    else
        if (strcmp(type, 'empiric') == 1)
            mp_func = @mp_empiric;
        elseif (strcmp(type, 'gauss') == 1)
            mp_func = @mp_gauss;
        elseif (strcmp(type, 'gaussi') == 1)
            mp_func = @mp_gaussi;
        elseif (strcmp(type, 'gammai') == 1)
            mp_func = @mp_gammai;
        else
            fprintf(2, ['\nValid Mutual Proximity type missing!\n'...
                'Use: Dmp = mutual_proximity(D, ''empiric''|'...
                '''gauss''|''gaussi''|''gammai'');\n\n']);
            Dmp = [];
            return;
        end
    end
    
    Dmp = mp_func(D);
end    


function Dmp = mp_empiric(D)

    n = size(D, 1);
    
    Dmp_list = cell(size(D, 1), 1);
    for i = 1:(n-1)
    
        % select only finite distances for MP
        j_idx = i+1:n;
        j_len = length(j_idx);
    
        dI = repmat(D(i, :), j_len, 1); 
        dJ = D(j_idx, :); 
        d = repmat(D(j_idx, i), 1, n); 
    
        sIJ_intersect = sum((dI > d) & (dJ > d), 2); 
    
        sIJ_overlap = 1 - (sIJ_intersect / n); 
        Dmp_list{i} = sIJ_overlap;
    end 
    
    Dmp = zeros(size(D), class(D));
    for i = 1:(n-1)
        j_idx = i+1:n;
        Dmp(i, j_idx) = Dmp_list{i}';
        Dmp(j_idx, i) = Dmp_list{i};
    end
end


function Dmp = mp_gaussi(D)

    mu = mean(D);
    sd = std(D);

    Dmp = zeros(size(D), class(D));
    n = length(D);
    
    for i=1:n
        
        j_idx = i+1:n;
        j_len = length(j_idx);
        
        p1 = 1 - local_normcdf(D(i, j_idx), ...
            repmat(mu(i), 1, j_len), repmat(sd(i), 1, j_len));
        p2 = 1 - local_normcdf(D(j_idx, i)', ...
            mu(j_idx), sd(j_idx));
        
        Dmp(i, j_idx) = 1 - p1.*p2;
        Dmp(j_idx, i) = Dmp(i, j_idx);

% Old Non-Vectorized Code (slow)
%         for j=i+1:n
%             
%             p1 = (1 - local_normcdf(D(i, j), mu(i), sd(i)));
%             p2 = (1 - local_normcdf(D(j, i), mu(j), sd(j)));
%             
%             Dmp(j, i) = 1 - p1*p2;
%             Dmp(i, j) = Dmp(j, i);
%                 
%         end

    end
end


function Dmp = mp_gauss(D)

    mu = mean(D);
    sd = std(D);
    
    % Ignore this warning for now
    warning off MATLAB:quadgk:MinStepSize
    epsmat = [100000*eps 0; 0 100000*eps];
    
    Dmp = zeros(size(D), class(D));
    n = size(D, 1);

    for i=1:n
        
        for j=i+1:n
            c = cov(D([i j], :)');
            x = [D(i, j) D(j, i)];
            m = [mu(i) mu(j)];
            
            p1 = local_normcdf(D(j, i), mu(i), sd(i));
            p2 = local_normcdf(D(j, i), mu(j), sd(j));
            try
                p12 = mvncdf(x, m, c);
            catch err
                if (strcmp(err.identifier,'stats:mvncdf:BadMatrixSigma'))
                    c = c + epsmat;
                    p12 = mvncdf(x, m, c);
                end
            end
            
            Dmp(j, i) = p1 + p2 - p12;
            Dmp(i, j) = Dmp(j, i);
                
        end
    end
end


function Dmp = mp_gammai(D)

    mu = mean(D);
    va = var(D);
    A = (mu.^2)./va;
    B = va ./ mu;
    
    Dmp = zeros(size(D), class(D));
    n = size(D, 1);
    
    for i=1:n
        
        j_idx = i+1:n;
        j_len = length(j_idx);
        
        p1 = 1 - local_gamcdf(D(i, j_idx), ...
            repmat(A(i), 1, j_len), repmat(B(i), 1, j_len));
        p2 = 1 - local_gamcdf(D(j_idx, i)', ...
            A(j_idx), B(j_idx));
        
        Dmp(i, j_idx) = 1 - p1.*p2;
        Dmp(j_idx, i) = Dmp(i, j_idx);        

% Old Non-Vectorized Code (slow)
%         for j=i+1:length(distm)
%             
%             a2 = (mu(j)^2)/va(j);
%             b2 = va(j)/mu(j);
%             
%             p1 = (1 - local_gamcdf(D(j, i), a1, b1));
%             p2 = (1 - local_gamcdf(D(j, i), a2, b2));
%             Dmp(j, i) = 1 - p1*p2;
%             
%             Dmp(i, j) = Dmp(j, i);
%                 
%         end
    end
end


function p = local_normcdf(x, mu, sd)
    z = (x-mu) ./ sd;
    p = 0.5 * erfc(-z ./ sqrt(2));
end


function p = local_gamcdf(x, a, b)
    a(a < 0) = NaN;
    b(b <= 0) = NaN;
    x(x < 0) = 0;
    z = x ./ b;
    p = gammainc(z, a);
end
