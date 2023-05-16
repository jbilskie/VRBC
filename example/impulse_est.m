function [h,t] = impulse_est(sync_data, preamble, t_preamble, VRBC, plot_yn, N)
    preamble_interp = interp1(t_preamble,preamble,(VRBC.PRI:VRBC.PRI:length(preamble)*VRBC.PRI));
    data = sync_data;
    % Just get preamble data
    y = (data(1:length(preamble_interp)));
    % Cross correlation
    rxy = xcorr(preamble_interp,y); 
    Ryx = fliplr(rxy(1:length(y))); 
    Ryx = Ryx(1:length(y));
    t = (VRBC.PRI:VRBC.PRI:length(Ryx)*VRBC.PRI);
    h = Ryx;
    h = 2.*h./length(h);
    h = movmean(h,N);
    % Plot
    if plot_yn == 1
        f = linspace(-VRBC.PRF/2,VRBC.PRF/2,1024);
        figure()
        subplot(2,1,1)
        hold on
        plot(t,h)
        title('Impulse Response Estimate')
        xlabel('Time [sec]'), ylabel('Amplitude')
        grid on
        subplot(2,1,2)
        hold on
        plot(f,10*log10(abs(fftshift(fft(h,1024)))))
        title('Frequency Response Estimate')
        xlabel('Frequency [Hz]'), ylabel('Magnitude [dB]')
        grid on
    end
end
