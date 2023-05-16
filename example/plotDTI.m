function plotDTI(cube, STFT, target)
    vWinType = @taylorwin;
    vWin = window(vWinType, STFT.window_len, 10, -50);
    % Look at a specific range bin
    data = select_range(cube.data, cube.config, cube.ranges, target.location, cube.NFFTrange);
    data = mean(data, 1);
    % Process data
    output_data = [];
    start_i = 1;
    end_i = STFT.window_len;
    complete = 0;
    while complete==0
        processed = (data(:,start_i:end_i));
        processed = permute(processed,[2,1]).*vWin;
        processed = fftshift(fft(processed, cube.NFFTvelocity, 1), 1);
        output_data = [output_data, processed];
        start_i = start_i+STFT.window_len-STFT.overlap;
        end_i = end_i+STFT.window_len-STFT.overlap;
        if end_i > length(data(1,:))
            complete = 1;
        end
    end
    % Plot
    mag = abs(output_data);
    norm_processed_out = mag./max(max(mag));
    intensities = 20*log10(norm_processed_out);
    time_bounds = [0, length(data(1,:))/cube.config.numChirps*cube.config.framePeriodicity*1e-3];
    doppler_bounds = [-max(cube.dopplers), max(cube.dopplers)];
    figure()
    imagesc(time_bounds, doppler_bounds, intensities(:,:))
    ylabel('Doppler [Hz]'), xlabel('Time [sec]')
    title('Dopper Time Intensity Plot')
    colormap jet
    caxis([-60,0]), colorbar
end

function processed = select_range(data, config, ranges, location, NFFTrange)
    rWinType = @taylorwin;
    rWin = window(rWinType, config.numADCSamples, 10, -50);
    % Process
    processed = rWin.*data;
    processed = fft(processed, NFFTrange, 1);
    diff = abs(ranges-location);
    diff_sort = sort(diff);
    index = find(diff==diff_sort(1));
    processed = processed(index,:);
end
