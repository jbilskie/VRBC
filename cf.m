function out = cf(in, num_chirp_for_cf, L)
    % Variable num_chirp_for_cf is the number of chirps averaged over for 
    % estimating the constant value leading up to the current chirp ---
    % This assumes 0 Hz clutter filtering (stationary clutter and target
    % platform.
    [m,lk] = size(in);
    K = ceil(lk/L);
    copy_in = in;
    out = zeros(size(in));
    for g = num_chirp_for_cf+1:lk
        curr_chirp = copy_in(:,g);
        prev_chirp_ave = mean(copy_in(:,g-num_chirp_for_cf:g-1),2);
        out(:,g) = curr_chirp-prev_chirp_ave;
    end
end
