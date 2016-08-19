function di = goodman_kruskal(D, classes)
% Computes the Goodman-Kruskal clustering index.
%
% This file is part of the HUB TOOLBOX available at
% http://ofai.at/research/impml/projects/hubology.html
% https://github.com/OFAI/hub-toolbox-matlab/
% (c) 2013, Dominik Schnitzer <dominik.schnitzer@ofai.at>
% (c) 2016, Roman Feldbauer <roman.feldbauer@ofai.at>
%
% The Goodman-Kruskal index is a
% clustering quality measure that relates the number of concordant ($Q_c$) and
% discordant (Q_d) tuples (d_{i,j}, d_{k,l}) of a distance matrix.
%  * A tuple is concordant if its items i, j are from the same class,
%    items k, l are from different classes and d_{i,j} < d_{k,l}.
%  * A tuple is discordant if its items i, j are from the same class,
%    items k, l are from different classes and d_{i,j} > d_{k,l}.
%  * A tuple is not counted if it is neither concordant nor discordant,
%    that is, if d_{i,j} = d_{k,l}.
%
% The Goodman-Kruskal Index ($I_{GK}$) is defined as:
% I_{GK} = \frac{Q_c - Q_d}{Q_c + Q_d}.
%
% I_{GK} is bounded to the interval [-1, 1], and the higher I_{GK}, the
% more concordant and fewer discordant quadruples are present in the data set.
% Thus a large index value indicates a good clustering (in terms of
% pairwise stability.
%
% Usage:
%   goodman_kruskal(D, classes) - Where D is an NxN distance matrix and
%      classes is a vector with the class labels as integers.

    Qc = 0;
    Qd = 0;
    cls = unique(classes);
    n = length(classes);
    
    % D_kl pairs in different classes
    D_other = D(triu(repmat(classes', n, 1) ~= repmat(classes, 1, n)))';
    for c = 1:length(cls)
        
        sel = classes == cls(c);
        if (sum(sel) > 1)
            
            selD = false(size(D));
            selD(sel, sel) = true;
            % D_ij pairs within same class
            D_self = D(triu(selD, 1))';
        else
            % skip if there is only one item per class
            continue;
        end
        % D_kl pairs in different classes (D_other) are computed once for all c
        D_full = [D_self D_other];
        
        self_size = length(D_self);
        other_size = length(D_other);
        
        [~, full_idx] = sort(D_full);
        
        % Calc number of quadruples with equal distance
        n_equidistant = 0;
        sdf = sort(D_full);
        equi_mask = boolean(zeros(size(sdf)));
        % Positions with repeated values
        equi_mask(2:end) = sdf(2:end) == sdf(end:-1:2);
        equi_dist = sdf(equi_mask);
        % How often does each value occur in self/other:
        uniq_equi_dist = unique(equi_dist);
        for i = 1:length(uniq_equi_dist)
            dist = uniq_equi_dist(i);
            equi_arg = find(D_full == dist);
            self_equi = sum(equi_arg < self_size);
            other_equi = length(equi_arg) - self_equi;
            % Number of dc that are actually equal
            n_equidistant = n_equidistant + self_equi * other_equi;
        end

        % Calc number of concordant quadruples
        cc = 0;
        ccsize = other_size;
        for i = 1:length(full_idx)
            if (full_idx(i) <= self_size)
                cc = cc + ccsize;
            else
                ccsize = ccsize - 1;
            end
        end
        
        % Calc number of discordant quadruples
        dc = self_size*other_size - cc - n_equidistant;
        
        Qc = Qc + cc;
        Qd = Qd + dc;
    end
    
    % Calc Goodman-Kruskal's gamma
    if (Qc+Qd == 0)
        di = 0;
    else
        di = (Qc-Qd)/(Qc+Qd);
    end
    
end