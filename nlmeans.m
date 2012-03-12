%%%%%%%%%%
% NL means presented in 
% http://hal.archives-ouvertes.fr/docs/00/27/11/41/PDF/061602r.pdf
% 
% Uses [9x9] similarity window, compares to window of [21x21]
% 
%%%%%%%%%%
function [NL_result] = nlmeans(I, h, a)
    sim = 9; % size of similarity window
    sim_w = (sim-1)/2; % size of sim on one side
    locality = 21; % 'neighborhood' patch. locality
    neigh_w = (locality-1)/2; 
    %R    OFFSET = w + sim_w;
    % pad an image so that we can handle edges
    I2 = padarray(I, [sim_w, sim_w], 'replicate'); 
    % assumed to be a square
    [n ~] = size(I); 
    % pixel location in original image
    % first extract patches = w x h x sim^2
    patches = zeros(n, n, sim^2);
    for x = 1:numel(I)
        [r, c] = ind2sub(size(I), x);
        window = I2(r:r+2*sim_w, c:c+2*sim_w); % offset by sim_w
        % if mod(x, 100)==0
        %     imshow(window);
        % end
        patches(r,c, :) = window(:);
    end

    % for all neighrbors, compute the weight to get teh weighted
    % average of all neighbors
    % gaussian weighted euclidean distance
    %http://www.econ.upf.edu/~michael/stanford/maeb4.pdf
    % || v2 - v1 ||^2_a = \sum j=length(v2) a(v2_j - v1_j)^2    
    weightedEuclideanDist = @(v2s, v1, a) ...
                            sum(a*bsxfun(@minus,v2s,v1).^2, 2);
    NL_result = zeros(size(I));
    for x = 1:numel(I)
        [r, c] = ind2sub(size(I), x);
        ind_r = [r - neigh_w: r+neigh_w];
        ind_c = [c - neigh_w : c+neigh_w];
        %take care of the boundary cases by reflection
        ind_r(ind_r < 1) = 2-ind_r(ind_r<1);
        ind_r(ind_r > n) = 2*n - ind_r(ind_r>n);
        ind_c(ind_c < 1) = 2-ind_c(ind_c<1);
        ind_c(ind_c > n) = 2*n - ind_c(ind_c>n);
        ngbhs = patches(ind_r, ind_c, :);
        ngbhs = reshape(ngbhs, [locality.^2, sim.^2]);
        v1 = squeeze(patches(r, c, :))';
        dists = weightedEuclideanDist(ngbhs, v1, a);
        origin = find(dists==0);
        rhs = exp(-dists/h.^2);
        Z = sum(rhs);
        weights = rhs./Z;
        neighbor_pixel_val = I(ind_r, ind_c);
        neighbor_pixel_val = neighbor_pixel_val(:);
        NL_result(r,c) = sum(neighbor_pixel_val.*weights);
    end
end

