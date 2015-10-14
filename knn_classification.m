function [acc, corr, cmat] = knn_classification(D, classes, k)
% Performs a k-nearest neighbor classification experiment. If there is a
% tie, the nearest neighbor determines the class
%
% This file is part of the HUB TOOLBOX available at
% http://ofai.at/research/impml/projects/hubology.html
% (c) 2013, Dominik Schnitzer <dominik.schnitzer@ofai.at>
%
% Usage:
%   [acc, corr, cmat] = knn_classification(D, classes, k) - Use the distance
%      matrix D (NxN) and the classes and perform a k-NN experiment. The
%      classification accuracy is returned in acc. corr is a raw vector of the
%      correctly classified items. cmat is the confusion matrix. 
	
    acc = zeros(length(k), 1);
    corr = zeros(length(D), length(k));
    
    n = length(D);
    
    cl = sort(unique(classes));
    cmat = zeros(length(cl));
    
    for i = 1:n
        classes(i) = find(cl == classes(i));
    end
        
	for i = 1:n
        
        seed_class = classes(i);
        
		row = D(i, :);
		row(i) = +Inf;
		[tmp, idx] = sort(row);
        
        for j = 1:length(k)
        
            nn_class = classes(idx(1:k(j)));
            cs = histc(nn_class, cl);
            id = find(cs == max(cs));

            % "tie": use nearest neighbor
            if (length(id) > 1)
                if (seed_class == nn_class) 
                    acc(j) = acc(j) + 1/n;
                    corr(i,j) = 1;
                end
                cmat(seed_class, nn_class) = cmat(seed_class, nn_class) + 1;

            % majority vote        
            else
                if (cl(id) == seed_class)
                    acc(j) = acc(j) + 1/n;
                    corr(i,j) = 1;
                end
                cmat(seed_class, cl(id)) = cmat(seed_class, cl(id)) + 1;
            end
            
        end
        
	end

end
