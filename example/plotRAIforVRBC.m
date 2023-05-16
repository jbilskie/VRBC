function [intensities] = plotRAIforVRBC(data, angs, config, NFFTrange, VRBC_freqs, NFFTdop, plott, fn, sbplt)
    % Get Doppler Indexes for Analysis
    PRI = (config.profileCfg.idleTime+config.profileCfg.rampEndTime)*1e-6;
    dopps = linspace(-1/(2*PRI),1/(2*PRI),NFFTdop+1);
    dopps = dopps(1:end-1);
    dopp_indexes = [];
    for d = 1:length(VRBC_freqs)
        temp = find(abs(dopps-VRBC_freqs(d))==min(abs(dopps-VRBC_freqs(d))));
        dopp_indexes = [dopp_indexes, temp];
        temp = find(abs(dopps+VRBC_freqs(d))==min(abs(dopps+VRBC_freqs(d))));
        dopp_indexes = [dopp_indexes, temp];
    end
    % Get VRBC Comparison Range-Angle Intensities
    yesdatacube = zeros(NFFTrange,length(angs),NFFTdop);
    tempdata = data;
    [M,N,LK] = size(tempdata);
    vWinType = @taylorwin;
    vWin = window(vWinType, LK, 10, -50);
    rWinType = @taylorwin;
    rWin = window(rWinType, M, 10, -50);
    for theta = angs
        W = beamformer(N, 0, tempdata, theta, config.profileCfg.startFreq*1e9);
        temp = zeros(M,N);
        for chirp_ind = 1:LK
            temp(:,chirp_ind) = squeeze(tempdata(:,:,chirp_ind))*W;
        end
        processed = rWin.*temp;
        processed = fft(processed, NFFTrange, 1);
        processed = vWin.*permute(processed,[2,1]);
        processed = fftshift(fft(processed, NFFTdop, 1),1);
        yesdatacube(:,angs==theta,:) = processed.';
    end
    yesmag = abs(yesdatacube);
    rangeAngle = abs(max(yesmag(:,:,dopp_indexes),[],3));
    norm_processed_out_r = rangeAngle./max(max(rangeAngle));
    intensities = 10*log10(norm_processed_out_r);
    % Bounds
    range_bounds = [0,(config.profileCfg.outSampleRate*1e3*config.c)/(2*config.profileCfg.rampSlope*1e12)];
    if plott == 1
        figure(fn)
        subplot(sbplt(1),sbplt(2),sbplt(3))
        imagesc(angs, range_bounds, intensities)
        ylabel('Range [m]'), xlabel('Angle [degree]')
        title('Range Angle Intensity Plot')
        colormap jet
        caxis([-20,0]), colorbar
        caxis([max(max(intensities))-30, max(max(intensities))])
        ylim([0, range_bounds(end)/2])
    end
end
