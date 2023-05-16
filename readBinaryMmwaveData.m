%% Reads ADC capture data from ISWR14XX using HSDC-pro into MATLAB

% Requires:
% 1) Binary file containing radar data
%    binFile = 'C:\User\jmb221\Document\test.bin'
% 2) Paramters expected to be in a configure structure
%    % Basic Configuration Parameters
%    config.numTX = 1;               % Number of TX antennas used
%    config.numRX = 1;               % Number of RX antennas used
%    config.numFrames = 160;         % Number of frames collected
%    config.numADCSamples = 530;     % Number of ADC samples per chirp
%    config.numChirps = 245;         % Number of chirps per frame
%    config.isreal = 0;              % Whether data is real or complex
%    config.numADCBits = 16;         % Number of ADC bits
%    config.c = 300000000;           % Speed of light
%    config.framePeriodicity = 16.1; % Frame period in [ms]
%    config.BW = 4000;               % Chirp bandwaidth [MHz]
%    % Chirp Profile Configuration Parameters
%    config.profileCfg.idleTime = 10;        % Idle time [usec]
%    config.profileCfg.rampEndTime = 52.5;   % Ramp end time [usec]
%    config.profileCfg.rampSlope = 80;       % Ramp slope [usec]
%    config.profileCfg.startFreq = 77;       % Start frequency [GHz]
%    config.profileCfg.ADCStartTime = 5;     % ADC start time [usec]
%    config.profileCfg.outSampleRate = 6000; % ADC sample rate [ksps]

% Dimension of output radar cube is:
% (ADC_SAMPLE_IDX) X (RX_IDX) X (CHIRPIDX) X (FRAME_IDX

% Note: Currently only setup for SIMO operation and interleaved data.

function data = readBinaryMmwaveData(binFile, config, radar, window_yn)
    fid = fopen(binFile, 'r');
    data = fread(fid, 'uint16');
    data = data-(data>=2^(config.numADCBits-1)).*2^config.numADCBits;
    fclose(fid);
    if config.isreal
        warning('Data is real. Please update code to handle.')
    else
        if radar == 77
            data = reshape(data ,config.numRX*2, []);
            data = (data(1:4,:)+1i*data(5:8,:));
        elseif radar == 60
            data = data(1:config.numRX*2*config.numADCSamples*config.numChirps*config.numFrames);
            data = reshape(data ,config.numRX*2*config.numADCSamples, []);
            rx0 = zeros(1,config.numFrames*config.numChirps*config.numADCSamples);
            rx1 = zeros(1,config.numFrames*config.numChirps*config.numADCSamples);
            rx2 = zeros(1,config.numFrames*config.numChirps*config.numADCSamples);
            rx3 = zeros(1,config.numFrames*config.numChirps*config.numADCSamples);
            for frame = 1:config.numFrames*config.numChirps
                for bit = 1:config.numADCSamples/2
                    rx0(((bit-1)*2+(1:2))+((frame-1)*config.numADCSamples)) = data((bit-1)*4+(1:2),frame)+1i*data((bit-1)*4+(3:4),frame);
                    rx1(((bit-1)*2+(1:2))+((frame-1)*config.numADCSamples)) = data((bit-1)*4+(1:2)+(2*config.numADCSamples),frame)+1i*data((bit-1)*4+(3:4)+(2*config.numADCSamples),frame);
                    rx2(((bit-1)*2+(1:2))+((frame-1)*config.numADCSamples)) = data((bit-1)*4+(1:2)+(4*config.numADCSamples),frame)+1i*data((bit-1)*4+(3:4)+(4*config.numADCSamples),frame);
                    rx3(((bit-1)*2+(1:2))+((frame-1)*config.numADCSamples)) = data((bit-1)*4+(1:2)+(6*config.numADCSamples),frame)+1i*data((bit-1)*4+(3:4)+(6*config.numADCSamples),frame);
                end
            end
            data = [rx0;rx1;rx2;rx3];
        else
            warning('Radar type is unrecognized. Please update code to handle.')
        end
    end
    data = data(:,1:config.numADCSamples*config.numChirps*config.numFrames);
    data = reshape(data, config.numRX, config.numADCSamples,config.numChirps, config.numFrames);
    data = permute(data, [2,1,3,4]);
    if window_yn == 1
        hann_win = hanning(size(data(:,1,1,1),1));
        hann_win = hann_win/sum(hann_win)/config.numChirps;
        for rx_i = 1:config.numRX
            for chirp_i = 1:config.numChirps
                for frame_i = 1:config.numFrames
                    data(:,rx_i,chirp_i,frame_i) = hann_win.*(data(:,rx_i,chirp_i,frame_i));
                end
            end
        end
    end
end
