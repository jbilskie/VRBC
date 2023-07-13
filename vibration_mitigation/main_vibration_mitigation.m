clc
clear
close all

%% system parameters
T = 1.80e-4;
beta = 4e9/T;
num_clt_pts = 1;
fs = 1e6;
full_ts = 1/fs:1/fs:0.05;
PRI = 1/5000;
f0 = 77e9;
lambda = 3e8/f0;
% Array Parameters
array_pos = -5*lambda/2:lambda/4:5*lambda/2;
% Radar Movement Parameters
radar_vel= [0,0];   % [velocity in m/sec, direction of movement]
% radar_vib = 1e-5.*lowpass(randn(size(full_ts)),180,fs)+...
%             5e-5.*sin(2*pi*60.*full_ts+2*pi*rand(1))+...
%             7e-5.*sin(2*pi*25.*full_ts+2*pi*rand(1));

seed1 = 0.3 % CONTROL

radar_vib = 1e-3.*sin(2*pi*50.*full_ts+2*pi*seed1);

% Clutter Point Parameters
var_clts = 10*ones(num_clt_pts,1);
range_clts = (fs*3e8)/(2*beta).*rand(num_clt_pts,1);

range_clts = 3.5 % CONTROL

az_angs_clts = 180*rand(num_clt_pts,1)-90;
% az_angs_clts = 120*rand(num_clt_pts,1)-60; % constrained
az_angs_clts = az_angs_clts(randperm(num_clt_pts));

az_angs_clts = 30 % CONTROL

el_angs_clts = zeros(length(az_angs_clts),1); % temporarily use zeros

%% Generate Simulated Datacubes of num_clt_pts stationary Clutter
cd 'C:\Users\James\Documents\Yr1PhDResearch\simulationCode\'
% recall simulate return phases are randomized
clt_out = simulateCltReturn(radar_vel, radar_vib, var_clts, range_clts, az_angs_clts, full_ts, beta, f0, PRI, T, array_pos);
cd 'C:\Users\James\Documents\Yr1PhDResearch\vibration_mitigation\'

% Add AWGN Noise
var_n = 1; % scale to signal amplitude
n_out = (sqrt(var_n)/2).*(randn(size(clt_out))+1j*randn(size(clt_out)));

% Full Datacube
out = clt_out + n_out;
out = clt_out;

% Radarcube Dimensions
[m,n,l] = size(out);

%% Range vs Angle for a given Set of Chirps
% Plot Bounds
dopps = [-0.5 0.5]*(1/PRI);
rangs = [0 (fs*3e8)/(2*beta)];
angs = [-90 90];
chirps = 1:100;
% Beamform
angs_span = linspace(angs(1),angs(2),100);
bf_out = zeros(length(angs_span),n,l);
for ang_i = 1:length(angs_span)
    cd 'C:\Users\James\Documents\Yr1PhDResearch\simulationCode'
    W = beamformer(array_pos, angs_span(ang_i), 3e8/f0);
    cd 'C:\Users\James\Documents\Yr1PhDResearch\vibration_mitigation\'
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

%% clutter phase estimation
phase = zeros(num_clt_pts,l);

for M = 1:num_clt_pts % typo in paper
    % beamform to clutter point M
    theta = az_angs_clts(M); % [degrees]
    bf_out = zeros(n,l);
    cd 'C:\Users\James\Documents\Yr1PhDResearch\simulationCode'
    W = beamformer(array_pos, theta, 3e8/f0);
    cd 'C:\Users\James\Documents\Yr1PhDResearch\vibration_mitigation\'
    for chirp_ind = 1:l
        bf_out(:,chirp_ind) = W.'*out(:,:,chirp_ind);
    end
    
    % Range FFT
    r_out = fft(bf_out(:,:),[],1);

    % Look at a specific range bin
    ranges = linspace(0,(fs*3e8)/(2*beta), n); % is this correct?
    diff = abs(ranges-range_clts(M));
    diff_sort = sort(diff);
    index = find(diff==diff_sort(1));

    s_vib = r_out(index,:); % this should be a complex scalar?

    % Frame Phase Estimation
    testimag = imag(s_vib);
    testreal = real(s_vib);
    test = imag(s_vib)./real(s_vib);
    phi_m = atan(imag(s_vib)./real(s_vib));
    phase(M,:) = phi_m;

end

%% estimation of global vibration effect
y_glob = zeros(num_clt_pts,l);
f_h = (2*radar_vel(1))/lambda;

for M = 1:num_clt_pts % typo in paper
    y_globM = (phase(M,:) - 2*pi*f_h*T) / (cosd(az_angs_clts(M))*cosd(el_angs_clts(M))); % is T "chirp duration"?
    y_glob(M,:) = y_globM;
end  
y_globavg = (1/num_clt_pts) * sum(y_glob,1);

%% vibration mitigation
% y_beam = [];
% s_beam = [];

% for b = 1:B
%     y_beam(b) = y_glob(b)*cos(theta_az(m))*cos(theta_el(m)); % project onto bth beam
%     s_beam(b) = 
% end

%% plots
% range-Angle of return

% Fig 4: Doppler spectrum with and without simulated vibration

% Fig 5: vibration phase, estimated phase, estimation error,
figure; plot(radar_vib);
figure; plot(y_globavg);

% Fig. 6: Doppler spectrum of no vibration and mitigated

% Fig. 7: radar platform vibration over time

