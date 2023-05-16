function [phase_data] = phaseEst(datacube, VRBC, range, num_chrps, nfft, plot_yn, clt)
    % Look at a specific range bin
    ranges = linspace(0,(VRBC.fs*3e8)/(2*VRBC.beta), nfft);
    rWinType = @taylorwin;
    rWin = window(rWinType, length(datacube(:,1)), 10, -50);
    processed = rWin.*datacube;
    processed = fft(processed, nfft, 1);
    diff = abs(ranges-range);
    diff_sort = sort(diff);
    index = find(diff==diff_sort(1));
    dataraw = processed(index,:);
    % Estimate Phase and Clutter Filter
    r = 0.95;
    b = [1, -2*cos(0), 1];
    a = [1, -2*r*cos(0), r^2];
    data = [];
    L = num_chrps;
    for k = 1:length(datacube(1,:))/L
        temp_data = unwrap(angle(dataraw(:,((k-1)*L)+1:k*L)));
        if clt == 1
            temp_data = filtfilt(b,a,temp_data);
        end
        data = [data, temp_data];
    end
    phase_data = data;
    if plot_yn == 1
        figure()
        plot(VRBC.PRI:VRBC.PRI:VRBC.PRI*length(datacube(1,:)), phase_data)
        ylabel('Amplitude'), xlabel('Time [sec]')
        title('Phase Signal')
    end
end
