function new_cube = interp_btwn_frames(datacube, config)
    [M,N,L,K] = size(datacube);
    
    datacube = reshape(datacube, M, N, L*K);
    
    PRI = (config.profileCfg.idleTime+config.profileCfg.rampEndTime)*1e-6;
    frameT = config.framePeriodicity*1e-3;
    
    t_slow = repmat(PRI:PRI:config.numChirps*PRI,1,config.numFrames)+...
             repelem((0:frameT:(config.numFrames-1)*frameT),config.numChirps);
    t_slow_full = PRI:PRI:(config.numFrames*config.framePeriodicity*1e-3);
    new_cube = zeros(M,N,length(t_slow_full));
    for n = 1:N
        for m = 1:M
            new_cube(m,n,:) = interp1(t_slow,permute(datacube(m,n,:),[1,3,2]),t_slow_full);
        end
    end
end
