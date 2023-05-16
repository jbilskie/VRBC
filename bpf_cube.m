function bpf_datacube = bpf_cube(bf_datacube, range, VRBC, width) 
    [m,lk] = size(bf_datacube);
    r_est_freq = VRBC.beta*2*range/3e8;
    r_est_freq = round(r_est_freq/VRBC.PRF)*VRBC.PRF;
    d1 = designfilt('bandpassfir', 'FilterOrder', floor(VRBC.samps_per_chirp/4), ...
                    'CutoffFrequency1', r_est_freq-(width*VRBC.PRF/2), ...
                    'CutoffFrequency2', r_est_freq+(width*VRBC.PRF/2), ...
                    'SampleRate', VRBC.fs);
    bpf_datacube = [];
    for chrp = 1:lk
        bpf_datacube = [bpf_datacube, filtfilt(d1,bf_datacube(:,chrp))]; 
    end
end