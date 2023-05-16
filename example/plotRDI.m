function plotRDI(data, config, plott, fn, sbplt)
    [M,LK] = size(data);
    vWinType = @taylorwin;
    vWin = window(vWinType, LK, 10, -50);
    rWinType = @taylorwin;
    rWin = window(rWinType, M, 10, -50);
    % Process data
    processed = rWin.*data;
    processed = fft(processed, M, 1);
    processed = vWin.*permute(processed,[2,1]);
    processed = fftshift(fft(processed, LK, 1),1);
    mag = abs(processed);
    intensities = 10*log10(mag);
    PRI = (config.profileCfg.idleTime+config.profileCfg.rampEndTime)*1e-6;
    range_bounds = [0,(config.profileCfg.outSampleRate*1e3*config.c)/(2*config.profileCfg.rampSlope*1e12)];
    dopp_bounds = [-1/(PRI*2), 1/(PRI*2)];
    if plott == 1
        figure(fn)
        subplot(sbplt(1),sbplt(2),sbplt(3))
        imagesc(dopp_bounds, range_bounds, intensities.')
        ylabel('Range [m]'), xlabel('Doppler [Hz]')
        title('Range Doppler Intensity Plot')
        colormap jet
        colorbar, caxis([max(max(intensities))-50, max(max(intensities))])
        ylim([0, range_bounds(end)/2])
    end
end
