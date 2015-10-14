function Dls = local_scaling(D, k, type)
% Transforms the given distance matrix into new one using local scaling [1]
% with the given neighborhood radius k. There are two types of local
% scaling methods implemented. The original one and NICDM, both reduce
% hubness in distance spaces, similarly to Mutual Proximity.
%
% This file is part of the HUB TOOLBOX available at
% http://ofai.at/research/impml/projects/hubology.html
% (c) 2013, Dominik Schnitzer <dominik.schnitzer@ofai.at>
%
% Usage:
%   Dls = local_scaling(D, k, type) - Applies local scaling to the distance
%      matrix D (NxN). The parameter k sets the neighborhood radius. type
%      the scaling type. The scaled distance matrix is returned.
%
% Possible types (type parameter):
%   'original': Original Local Scaling using the distance of the k'th
%      nearest neighbor.
%   'nicdm': Local Scaling using the average distance of the k nearest
%      neighbors.
%
%
% [1] Local and global scaling reduce hubs in space, 
% Schnitzer, Flexer, Schedl, Widmer, Journal of Machine Learning Research 2012

    if (nargin < 3)
        ls_func = @ls_k;
        fprintf('No Local Scaling type given. Using: ''original''\n');
        fprintf('For more stable results use: ''nicdm''\n');
    end
    if (nargin < 2)
        k = 7;
        fprintf('No neighborhood radius given. Using k=7\n');
    end
    if (nargin >= 3)
        if (strcmp(type, 'original') == 1)
            ls_func = @ls_k;
        elseif (strcmp(type, 'nicdm') == 1)
            ls_func = @ls_nicdm;
        else
            fprintf(2, ['\nValid Local Scaling type missing!\n'...
                'Use: Dls = local_scaling(D, ''original''|''nicdm'');\n\n']);
            Dls = [];
            return;        
        end
    end
    
    Dls = ls_func(D, k);

end


function Dls = ls_k(D, k)

    r = zeros(length(D), 1);
    for i = 1:length(D)
        di = D(i, :);
        di(i) = +Inf;
        [tmp, nn] = sort(di);
        r(i) = di(nn(k));
    end
    
    n = size(D, 1);
    Dls = zeros(size(D), class(D));
    for i = 1:n
        for j = i+1:n
            Dls(i, j) = D(i, j) / sqrt( r(i) * r(j) );
            Dls(j, i) = Dls(i, j);
        end
    end
end


function Dnicdm = ls_nicdm(D, k)

    r = zeros(length(D), 1);
    for i = 1:length(D)
        di = D(i, :);
        di(i) = +Inf;
        [tmp, nn] = sort(di);
        r(i) = mean(di(nn(1:k)));
    end
    rg = local_geomean(r);
    
    Dnicdm = zeros(size(D), class(D));
    for i = 1:length(D)
        for j = i+1:length(D)
            Dnicdm(i, j) = (rg*D(i, j)) / sqrt( r(i) * r(j));
            Dnicdm(j, i) = Dnicdm(i, j);
        end
    end
end

function m = local_geomean(x)
    m = exp(sum(log(x))./length(x));
end
