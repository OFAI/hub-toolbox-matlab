function hubness_analysis(D, classes, vectors)
% Performs a quick hubness analysis with all the functions provided in this
% toolbox.
%
% This file is part of the HUB TOOLBOX available at
% http://ofai.at/research/impml/projects/hubology.html
% https://github.com/OFAI/hub-toolbox-matlab/
% (c) 2013, Dominik Schnitzer <dominik.schnitzer@ofai.at>
% (c) 2016, Roman Feldbauer <roman.feldbauer@ofai.at>
%
% Usage:
%  hubness_analysis() - Loads the example data set and performs the
%     analysis
%
%  hubness_analysis(D, classes, vectors) - Uses the distance matrix D (NxN)
%     together with an optional class labels vector (classes) and the
%     original (optional) data vectors (vectors, Pts x Dims) to perform a
%     full hubness analysis

    haveClasses = false;
    haveVectors = false;
    if (nargin == 0)
        [D, classes, vectors] = load_dexter();
        haveClasses = true;
        haveVectors = true;
    elseif (nargin == 1)
        % all ok, just analyze D
    elseif (nargin == 2)
        haveClasses = true;
    else
        haveClasses = true;
        haveVectors = true;
    end
        

    n = size(D,1);
    [Sn5, tmp, Nk5] = hubness(D, 5);
    fprintf('\nHubness Analysis\n\n');
    
    fprintf('ORIGINAL DATA:\n');
    fprintf('data set hubness (S^n=5)                : %.2f\n', Sn5);
    fprintf('%% of anti-hubs at k=5                   : %.2f%%\n',...
        100*sum(Nk5==0)/n);
    fprintf('%% of k=5-NN lists the largest hub occurs: %.2f%%\n',...
        100*max(Nk5)/n);
    if (haveClasses == true)
    fprintf('k=5-NN classification accuracy          : %.2f%%\n',...
        100*knn_classification(D, classes, 5));
    fprintf('Goodman-Kruskal index (higher=better)   : %.3f\n',...
        goodman_kruskal(D, classes));
    else
    fprintf('k=5-NN classification accuracy          : No classes given\n');
    fprintf('Goodman-Kruskal index (higher=better)   : No classes given\n');
    end
    if (haveVectors == true)
    fprintf('original dimensionality                 : %d\n', size(vectors, 2));
    fprintf('intrinsic dimensionality estimate       : %d\n',...
        round(intrinsic_dim(vectors)));
    else
    fprintf('original dimensionality                 : No vectors given\n');
    fprintf('intrinsic dimensionality estimate       : No vectors given\n');
    end  
    
    
    fprintf('\nMUTUAL PROXIMITY (Empiric/Slow):\n');
    Dn = mutual_proximity(D, 'empiric');
    [Sn5, tmp, Nk5] = hubness(Dn, 5);
    fprintf('data set hubness (S^n=5)                : %.2f\n', Sn5);
    fprintf('%% of anti-hubs at k=5                   : %.2f%%\n',...
        100*sum(Nk5==0)/n);
    fprintf('%% of k=5-NN lists the largest hub occurs: %.2f%%\n',...
        100*max(Nk5)/n);
    if (haveClasses == true)
    fprintf('k=5-NN classification accuracy          : %.2f%%\n',...
        100*knn_classification(Dn, classes, 5));
    fprintf('Goodman-Kruskal index (higher=better)   : %.3f\n',...
        goodman_kruskal(Dn, classes));
    else
    fprintf('k=5-NN classification accuracy          : No classes given\n');
    fprintf('Goodman-Kruskal index (higher=better)   : No classes given\n');
    end

%      fprintf('\nMUTUAL PROXIMITY (Gauss):\n');
%     Dn = mutual_proximity(D, 'gauss');
%     [Sn5, tmp, Nk5] = hubness(Dn, 5);
%     fprintf('data set hubness (S^n=5)                : %.2f\n', Sn5);
%     fprintf('%% of anti-hubs at k=5                   : %.2f%%\n',...
%         100*sum(Nk5==0)/n);
%     fprintf('%% of k=5-NN lists the largest hub occurs: %.2f%%\n',...
%         100*max(Nk5)/n);
%     if (haveClasses == true)
%     fprintf('k=5-NN classification accuracy          : %.2f%%\n',...
%         100*knn_classification(Dn, classes, 5));
%     fprintf('Goodman-Kruskal index (higher=better)   : %.3f\n',...
%         goodman_kruskal(Dn, classes));
%     else
%     fprintf('k=5-NN classification accuracy          : No classes given\n');
%     fprintf('Goodman-Kruskal index (higher=better)   : No classes given\n');
%    end
% 
%     fprintf('\nMUTUAL PROXIMITY (Gaussi):\n');
%     Dn = mutual_proximity(D, 'gaussi');
%     [Sn5, tmp, Nk5] = hubness(Dn, 5);
%     fprintf('data set hubness (S^n=5)                : %.2f\n', Sn5);
%     fprintf('%% of anti-hubs at k=5                   : %.2f%%\n',...
%         100*sum(Nk5==0)/n);
%     fprintf('%% of k=5-NN lists the largest hub occurs: %.2f%%\n',...
%         100*max(Nk5)/n);
%     if (haveClasses == true)
%     fprintf('k=5-NN classification accuracy          : %.2f%%\n',...
%         100*knn_classification(Dn, classes, 5));
%     fprintf('Goodman-Kruskal index (higher=better)   : %.3f\n',...
%         goodman_kruskal(Dn, classes));
%     else
%     fprintf('k=5-NN classification accuracy          : No classes given\n');
%     fprintf('Goodman-Kruskal index (higher=better)   : No classes given\n');
%     end

    fprintf('\nMUTUAL PROXIMITY (Gammai):\n');
    Dn = mutual_proximity(D, 'gammai');
    [Sn5, tmp, Nk5] = hubness(Dn, 5);
    fprintf('data set hubness (S^n=5)                : %.2f\n', Sn5);
    fprintf('%% of anti-hubs at k=5                   : %.2f%%\n',...
        100*sum(Nk5==0)/n);
    fprintf('%% of k=5-NN lists the largest hub occurs: %.2f%%\n',...
        100*max(Nk5)/n);
    if (haveClasses == true)
    fprintf('k=5-NN classification accuracy          : %.2f%%\n',...
        100*knn_classification(Dn, classes, 5));
    fprintf('Goodman-Kruskal index (higher=better)   : %.3f\n',...
        goodman_kruskal(Dn, classes));
    else
    fprintf('k=5-NN classification accuracy          : No classes given\n');
    fprintf('Goodman-Kruskal index (higher=better)   : No classes given\n');
    end

    
    fprintf('\nLOCAL SCALING (Original, k=10):\n');
    Dn = local_scaling(D, 10, 'original');
    [Sn5, tmp, Nk5] = hubness(Dn, 5);
    fprintf('data set hubness (S^n=5)                : %.2f\n', Sn5);
    fprintf('%% of anti-hubs at k=5                   : %.2f%%\n',...
        100*sum(Nk5==0)/n);
    fprintf('%% of k=5-NN lists the largest hub occurs: %.2f%%\n',...
        100*max(Nk5)/n);
    if (haveClasses == true)
    fprintf('k=5-NN classification accuracy          : %.2f%%\n',...
        100*knn_classification(Dn, classes, 5));
    fprintf('Goodman-Kruskal index (higher=better)   : %.3f\n',...
        goodman_kruskal(Dn, classes));
    else
    fprintf('k=5-NN classification accuracy          : No classes given\n');
    fprintf('Goodman-Kruskal index (higher=better)   : No classes given\n');
    end

    
    fprintf('\nSHARED NEAREST NEIGHBORS (k=10):\n');
    Dn = shared_nn(D, 10);
    [Sn5, tmp, Nk5] = hubness(Dn, 5);
    fprintf('data set hubness (S^n=5)                : %.2f\n', Sn5);
    fprintf('%% of anti-hubs at k=5                   : %.2f%%\n',...
        100*sum(Nk5==0)/n);
    fprintf('%% of k=5-NN lists the largest hub occurs: %.2f%%\n',...
        100*max(Nk5)/n);
    if (haveClasses == true)
    fprintf('k=5-NN classification accuracy          : %.2f%%\n',...
        100*knn_classification(Dn, classes, 5));
    fprintf('Goodman-Kruskal index (higher=better)   : %.3f\n',...
        goodman_kruskal(Dn, classes));
    else
    fprintf('k=5-NN classification accuracy          : No classes given\n');
    fprintf('Goodman-Kruskal index (higher=better)   : No classes given\n');
    end

    fprintf('\n');
end

function [D, classes, vectors] = load_dexter()

    fprintf('\nNO PARAMETERS GIVEN! Loading & evaluating DEXTER data set.\n\n');
    fprintf('DEXTER is a text classification problem in a bag-of-word\n');
    fprintf('representation. This is a two-class classification problem\n');
    fprintf('with sparse continuous input variables.\n');
    fprintf('This dataset is one of five datasets of the NIPS 2003 feature\n');
    fprintf('selection challenge.\n\n');

    fprintf('http://archive.ics.uci.edu/ml/datasets/Dexter\n\n');

    n = 300;
    dim = 20000;
        
    vectors = zeros(n, dim);
    classes = zeros(n, 1);
    
    
    fid = fopen('example_datasets/dexter_train.data', 'r');
    for i=1:n
        tline = fgetl(fid);
        d = sscanf(tline, '%d:%d');
        vectors(i,d(1:2:end)) = d(2:2:end);
    end
    fclose(fid);

    fid = fopen('example_datasets/dexter_train.labels', 'r');
    for i=1:n
        tline = fgetl(fid);
        d = sscanf(tline, '%d');
        classes(i) = d(1);
    end
    fclose(fid);
    
    D = cosine_distance(vectors);

end

function D = cosine_distance(x)
    xn = sqrt(sum(x.^2, 2));
    x = x ./ repmat(xn, 1, size(x, 2));
    D = 1 - x*x';
    D(D<0) = 0;
    D = triu(D, 1) + triu(D, 1)';
end

