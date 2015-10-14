function di = goodman_kruskal(D, classes)
% Computes the Goodman-Kruskal clustering index.
%
% This file is part of the HUB TOOLBOX available at
% http://ofai.at/research/impml/projects/hubology.html
% (c) 2013, Dominik Schnitzer <dominik.schnitzer@ofai.at>
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
    for c = 1:length(cls)
        
        sel = classes == cls(c);
        if (sum(sel) > 1)
            
            selD = false(size(D));
            selD(sel, sel) = true;
            D_self = D(triu(selD, 1))';
        else
            % skip if there is only one item per class
            continue;
        end
        D_other = D(sel, ~sel);
        D_other = D_other(:)';
        D_full = [D_self D_other];
        
        self_size = length(D_self);
        other_size = length(D_other);
        
        [tmp, full_idx] = sort(D_full);
        
        
        cc = 0;
        ccsize = other_size;
        for i = 1:length(full_idx)
            if (full_idx(i) <= self_size)
                cc = cc + ccsize;
            else
                ccsize = ccsize - 1;
            end
        end
        
        dc = self_size*other_size - cc;
        
        Qc = Qc + cc;
        Qd = Qd + dc;
    end
    
    if (Qc+Qd == 0)
        di = 0;
    else
        di = (Qc-Qd)/(Qc+Qd);
    end
    
end
