clear all, close all, clc

% Waveform Parameters
T = 1.80e-4;
beta = 4e9/T;
num_clt_pts = 20;
fs = 1e6;
full_ts = 1/fs:1/fs:0.05;
PRI = 1/5000;
f0 = 77e9;
lambda = 3e8/f0;
% Array Parameters
array_pos = -5*lambda/2:lambda/4:5*lambda/2;
% Radar Movement Parameters
radar_vel= [0,0];   % [velocity in m/sec, direction of movement]
radar_vib = 1e-5.*lowpass(randn(size(full_ts)),180,fs)+...
            5e-5.*sin(2*pi*60.*full_ts+2*pi*rand(1))+...
            7e-5.*sin(2*pi*25.*full_ts+2*pi*rand(1));
% Clutter Point Parameters
var_clts = 10*ones(num_clt_pts,1);
range_clts = (fs*3e8)/(2*beta).*rand(num_clt_pts,1);
az_angs_clts = 180*rand(num_clt_pts,1)-90;
az_angs_clts = az_angs_clts(randperm(num_clt_pts));
% Transponder Parameters
var_tar = 10;
range_tar = 4;
az_ang_tar = 0;
tar_vib = 7e-5.*sin(2*pi*250.*full_ts);

% Generate Simulated Datacubes of Clutter and VRBC Transponder Scatterer
clt_out = simulateCltReturn(radar_vel, radar_vib, var_clts, range_clts, az_angs_clts, full_ts, beta, f0, PRI, T, array_pos);
tar_out = simulateVRBCReturn(radar_vel, radar_vib, var_tar, range_tar, az_ang_tar, tar_vib, full_ts, beta, f0, PRI, T, array_pos);

% Add AWGN Noise
var_n = 10;
n_out = (sqrt(var_n)/2).*(randn(size(clt_out))+1j*randn(size(clt_out)));

% Plot Bounds
dopps = [-0.5 0.5]*(1/PRI);
rangs = [0 (fs*3e8)/(2*beta)];
angs = [-90 90];

% Full Datacube
out = tar_out + clt_out + n_out;

% Radarcube Dimensions
[m,n,l] = size(out);

%% Range vs Dopper at a given Angle
theta = az_ang_tar; % [degrees]
% Beamform 
bf_out = zeros(n,l);
W = beamformer(array_pos, theta, 3e8/f0);
for chirp_ind = 1:l
    bf_out(:,chirp_ind) = W.'*out(:,:,chirp_ind);
end
% Range FFT
r_out = fft(bf_out(:,:),[],1);
% Doppler FFT
d_out = fftshift(fft(r_out(:,:),[],2),2);
% Plot
figure()
imagesc(dopps, rangs, 10*log10(abs(d_out)))
title('Range vs Doppler')
xlabel('Doppler [Hz]'), ylabel('Range [m]')

%% Range vs Angle for a given Set of Chirps
chirps = 1:100;
% Beamform
angs_span = linspace(angs(1),angs(2),100);
bf_out = zeros(length(angs_span),n,l);
for ang_i = 1:length(angs_span)
    W = beamformer(array_pos, angs_span(ang_i), 3e8/f0);
    for chirp_ind = 1:l
        bf_out(ang_i,:,chirp_ind) = W.'*out(:,:,chirp_ind);
    end
end
% Range FFT
r_out = fft(bf_out(:,:,chirps),[],2);
% Capture all the Power at a Given Range/Angle
d_out = sum(abs(r_out),3); 
% Plot
figure()
imagesc(angs, rangs, 10*log10(d_out).')
title('Range vs Angle')
xlabel('Angle [degrees]'), ylabel('Range [m]')