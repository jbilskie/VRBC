clc
clear all 
close all

%% Setup and Pulling in Data
% 1) Binary file containing radar data
   binFile = 'test.bin';
% 2) Paramters expected to be in a configure structure
   % Basic Configuration Parameters
   config.numTX = 1;               % Number of TX antennas used
   config.numRX = 4;               % Number of RX antennas used
   config.numFrames = 200;         % Number of frames collected
   config.numADCSamples = 539;     % Number of ADC samples per chirp
   config.numChirps = 124;         % Number of chirps per frame
   config.isreal = 0;              % Whether data is real or complex
   config.numADCBits = 16;         % Number of ADC bits
   config.c = 3e8;                 % Speed of light
   config.framePeriodicity = 25;   % Frame period in [ms]
   config.BW = 3999.5;             % Chirp bandwidth [MHz]
   % Chirp Profile Configuration Parameters
   config.profileCfg.idleTime = 10;        % Idle time [usec]
   config.profileCfg.rampEndTime = 190;    % Ramp end time [usec]
   config.profileCfg.rampSlope = config.BW/config.profileCfg.rampEndTime;
   config.profileCfg.startFreq = 77;       % Start frequency [GHz]
   config.profileCfg.ADCStartTime = 10;    % ADC start time [usec]
   config.profileCfg.outSampleRate = 3000; % ADC sample rate [ksps]
   % Resultant VRBC Parameters
   VRBC.PRI = (config.profileCfg.idleTime+config.profileCfg.rampEndTime)*1e-6;
   VRBC.PRF = 1/VRBC.PRI;
   VRBC.lambda = config.c/(config.profileCfg.startFreq*1e9);
   VRBC.fs = config.profileCfg.outSampleRate*1e3;
   VRBC.f0 = config.profileCfg.startFreq*1e9;
   VRBC.beta = config.profileCfg.rampSlope*1e12;
   VRBC.T = VRBC.PRI-(config.profileCfg.idleTime*1e-6);
   VRBC.samps_per_chirp = round(VRBC.T*VRBC.fs);
   VRBC.samps_per_PRI = round(VRBC.PRI*VRBC.fs);
   VRBC.time_obs = config.numFrames*config.framePeriodicity*1e-3;
   VRBC.L = round((config.framePeriodicity*1e-3)/VRBC.PRI);
   VRBC.K = config.numFrames; 
   % Chosen VRBC Parameters
   VRBC.Tsym = 0.02;                    % Length of a symbol [sec]
   VRBC.M = 4;                          % Number of symbols [#]
   VRBC.sym_freqs = [150,350,450,750];  % FSK frequencies [Hz]
   VRBC.Lsym = VRBC.Tsym*VRBC.fs;       % Samples per symbol [#]
   VRBC.num_sym = (VRBC.time_obs-2)/VRBC.Tsym;
   VRBC.chirps_per_sym = VRBC.Tsym/VRBC.PRI;
   % Calculate Excitation Signals
   x = [];
   for m=1:VRBC.M
       t1 = 0:1/VRBC.fs:VRBC.Tsym;
       x = [x; sin(2*pi*VRBC.sym_freqs(m).*t1)];
   end
   VRBC.x = x;
   clear x t1 t2

% Read Data
[filepath, name, ext] = fileparts(binFile);
axes = axesParameters(config);
datacube = readBinaryMmwaveData(binFile, config, 77, 0);
% Desired cube dimensions
M = config.numADCSamples;
N = config.numRX;
L = config.numChirps;
K = config.numFrames;
% Actual cube dimensions
[m,n,l,k] = size(datacube);
if sum([M,N,L,K]~=[m,n,l,k]) ~= 0
    warning('Size of datacube is unexpected.')
end
% Interpolate
datacube = interp_btwn_frames(datacube, config);
datacube(isnan(datacube)) = 0;
datacube = reshape(datacube,m,n,[],k);
clear M N K L m n l k

[num_samps, num_elems, num_chrps, num_frames] = size(datacube);
[ranges, dopplers, velocities] = computeAxes(axes, num_samps, num_chrps*num_frames);

%% Target Detection Via Doppler-Constrained Detection in Range-Angle
% Plot Range vs Angle Intensity (RAI) Plots
angs = -90:5:90;
yes_times = 2.4:.1:2.5; % specify frame start times to plot RAI
frame_times = (config.framePeriodicity*1e-3).*(1:num_frames);
for f = 1:length(yes_times)
    yes_frame_ind = find(abs(frame_times-yes_times(f))==min(abs(frame_times-yes_times(f))));
    temp_cube = squeeze(datacube(:,:,:,yes_frame_ind));
    [intensities] = plotRAIforVRBC(temp_cube, angs, config, num_samps, VRBC.sym_freqs, num_chrps, 1, 1, [1,2,f]);
    title(sprintf('Range Angle Intensity Plot\nTime = %.3f-%.3f sec',frame_times(yes_frame_ind-1),frame_times(yes_frame_ind)))
end

% Provide/find angle (theta) and range of a detected transponder.
[range_i,theta_i] = find(intensities==max(max(intensities)));
theta = angs(theta_i);    % [degrees]
range = ranges(range_i);  % [meters]

%% Target Isolation with Beamforming and Range Filtering
new_datacube = reshape(datacube, num_samps, num_elems, num_chrps*num_frames);

% Beamform
bf_datacube = zeros(num_samps,num_chrps*num_frames);
W = beamformer(num_elems, axes, datacube, theta, config.profileCfg.startFreq*1e9);
for chirp_ind = 1:num_chrps*num_frames
    bf_datacube(:,chirp_ind) = new_datacube(:,:,chirp_ind)*W;
end
plotRDI(bf_datacube(:,(yes_frame_ind-1)*num_chrps+1:yes_frame_ind*num_chrps), config, 1, 2, [1,3,1])
title(sprintf('Range Doppler Intensity Plot\nAfter Beamforming'))
cl = caxis;

% Bandpass Filter
bpf_width = 1;
bpf_datacube = bpf_cube(bf_datacube, range, VRBC, bpf_width);
plotRDI(bpf_datacube(:,(yes_frame_ind-1)*num_chrps+1:yes_frame_ind*num_chrps), config, 1, 2, [1,3,2])
title(sprintf('Range Doppler Intensity Plot\nAfter Range Bandpass Filtering'))
caxis(cl)

% Clutter Filter
cf_datacube = cf(bpf_datacube, 10, num_chrps);
plotRDI(cf_datacube(:,(yes_frame_ind-1)*num_chrps+1:yes_frame_ind*num_chrps), config, 1, 2, [1,3,3])
title(sprintf('Range Doppler Intensity Plot\nAfter Clutter Filtering'))
caxis(cl)

%% Phase Estimation and Interpolation
% Frame Phase Estimation
phase_signal = phaseEst(bpf_datacube, VRBC, range, num_chrps, 539, 1, 0);
title('Extracted Phase Signal for Sync'), grid on
t_slow_full = VRBC.PRI:VRBC.PRI:(num_frames*config.framePeriodicity*1e-3);

%% Synchronization throguh Magnitude Square Coherence
t_preamble = VRBC.PRI:VRBC.PRI:1.5; % time vector for the known preamble
preamble = chirp(t_preamble,1,1.5,VRBC.PRF/2); % known preamble

% Magnitude Squared Coherence to Synchronize
[strt_time, ~, output_data, f_msc, temp_slice] = MST_Coherence(phase_signal, preamble, VRBC.PRF, length(preamble)-200, [0,VRBC.PRF/2]);
figure(), imagesc([0, t_slow_full(end)-1], f_msc, abs(output_data))
title('Magnitude Squared Coherence with Preamble Across Time')
xlabel('Time [sec]'), ylabel('Frequency [Hz]')
look_f = VRBC.sym_freqs;
figure()
for s = 1:length(look_f)
    ind = find(abs(f_msc-look_f(s))==min(abs(f_msc-look_f(s))));
    hold on,
    plot(linspace(0, t_slow_full(end)-1,length(output_data(1,:))),output_data(ind,:),'DisplayName',sprintf("Freq = %i Hz",round(f_msc(ind))))
end
title('Slices of MSC Frequencies Across Time')
xlabel('Time [sec]'), ylabel('Amplitude'), grid on
legend()

strt_slowtime_ind = find(min(abs(t_slow_full-strt_time)) == abs(t_slow_full-strt_time) );
sync_phase = phase_signal(strt_slowtime_ind:end);

%% Impulse Respnse Estimation
% Take Away DC Component
sync_phase = sync_phase - mean(sync_phase(end-1000:end));

% Impulse Response Estimation
[h,h_t] = impulse_est(sync_phase(1:length(preamble)), preamble, t_preamble, VRBC, 1, 1);

% Estimate Sync Data
f = VRBC.sym_freqs;
seq = reshape(dec2base('Hello World! I am working.',VRBC.M).',1,[]);
binseq = [];
for n = 1:length(seq)
    binseq = [binseq, str2num(seq(n))];
end
Tsym = VRBC.Tsym;
t = 1/VRBC.PRF:1/VRBC.PRF:1;
t2 = 1/VRBC.PRF:1/VRBC.PRF:Tsym;
z = [];
for n = 1:length(seq)
    z = [z, sin(2*pi*f(binseq(n)+1).*t2)];
end
signal = [preamble,zeros(1,0.5*VRBC.PRF),z];
est_sig = conv(signal,h);

% Plot the Observed Phase Estimate vs the Predicted Using our h(t) Estimate
figure()
plot(VRBC.PRI.*(1:length(est_sig)),est_sig),hold on,plot(VRBC.PRI.*(1:length(sync_phase)),sync_phase)
legend('Estimated Slow-Time Phase','Actual Slow-Time Phase')
title('Comparison of Slow-Time Phase Signals - No CF')
xlabel('Time [sec]'),ylabel('Amplitude')

%% Displacement Modeling
M_colors = ['r';'b';'g';'m'];
d = [];
for m = 1:VRBC.M
   new_x = interp1(0:1/VRBC.fs:VRBC.Tsym,VRBC.x(m,:),1/VRBC.PRF:1/VRBC.PRF:VRBC.Tsym);
   temp = [0,conv(new_x, h)];
   temp = interp1(1/VRBC.PRF.*(0:1:length(temp)-1),temp,1/VRBC.fs:1/VRBC.fs:VRBC.PRI*(length(temp)-1));
   d = [d, temp.'];
end
figure()
for m = 1:VRBC.M
   subplot(2,1,1)
   hold on, plot((1/VRBC.fs).*(1:length(temp)),d(:,m),M_colors(m))
end
xlabel('Time [sec]'), ylabel('Amplitude'), grid on
title('Symbol Displacements')

% Calculate how long it takes each symbol to decay out.
min_dis_ratio = 0.5;                        % minimum displacement ratio 
min_dis = min_dis_ratio*max(max(abs(d)));   % minimum displacement allowed
set_N = zeros(1,VRBC.M);
for m = 1:VRBC.M
    temp_i = find(abs(d(:,m))>min_dis);
    if isempty(temp_i)
        temp_i = VRBC.Lsym;
    end
    set_N(m) = ceil((temp_i(end)-VRBC.Lsym)/VRBC.Lsym);
end

% Possible sequences displacements.
N = max(set_N);
combs = [];
for n = 1:N
    combs_temp = dec2base(0:VRBC.M^n-1,VRBC.M)' - '0';
    combs_temp = combs_temp+1;
    combs_temp = [zeros(N-n,length(combs_temp(1,:)));combs_temp];
    combs = [combs, combs_temp];
end
combs = flipud(combs);

% Eliminate unnecessary contributions and delete repeated sequences.
for s = 1:length(combs(1,:))
    for bit = 1:length(combs(:,1))
        set_ignore = 1:VRBC.M;
        ign_ind = find(set_N<bit);
        set_ignore= set_ignore(ign_ind);
        if ismember(combs(bit,s),set_ignore)==1
            combs(bit,s) = 0;
        end
    end
end
combs = [zeros(N,1),combs];
combs = unique(combs.','rows');
states = fliplr(combs);
S = length(states(:,1));
M_sym = ones(S,VRBC.M);

% State Transition Matrix
A = zeros(S);
for s = 1:S
    for m = 1:VRBC.M
        temp_seq = [states(s,2:end),m];
        for bit = 1:N
            set_ignore = 1:VRBC.M;
            ign_ind = find(set_N<bit);
            set_ignore= set_ignore(ign_ind);
            if ismember(temp_seq(end-bit+1),set_ignore)==1
                temp_seq(end-bit+1) = 0;
            end
        end
        [~, ind]=ismember(temp_seq,states,'rows');
        A(s,ind) = A(s,ind) + (1/VRBC.M);
    end
end

% State probabilities at each symbol interval
num_sym = length(binseq);
pis = zeros(S,num_sym+2);
pis(:,1) = [1; zeros(S-1,1)];
for s = 2:num_sym+1
    pis(:,s) = pis(:,s-1).'*A;
end

% Calculate all displacements
d_state_est = zeros(VRBC.Lsym,nnz(M_sym));
state_ind_vec = ceil(find(M_sym.'==1)/VRBC.M);
ending_sym_vec = repmat(1:VRBC.M,1,S);
ending_sym_vec = ending_sym_vec(M_sym.'==1);
for ss = 1:nnz(M_sym)
    s = state_ind_vec(ss);
    for p = 1:length(states(1,:))
        if states(s,end-p+1) ~= 0
            d_state_est(:,ss) = d_state_est(:,ss)+ d((p)*VRBC.Lsym+(1:VRBC.Lsym),states(s,end-p+1));
        end
    end
    d_state_est(:,ss) = d_state_est(:,ss) + d(1:VRBC.Lsym,ending_sym_vec(ss));
end
for ss = 1:nnz(M_sym)
    subplot(2,VRBC.M,VRBC.M+ending_sym_vec(ss))
    hold on, plot((1/VRBC.fs).*(1:length(d_state_est(:,ss))),d_state_est(:,ss),M_colors(ending_sym_vec(ss)))
end    
for m = 1:VRBC.M
    subplot(2,VRBC.M,VRBC.M+m)
    xlabel('Time [sec]'), ylabel('Amplitude'), grid on
    title(sprintf('State Disp. (ending in %i)',m))
end

% Calculate delays.
tau_tar_est = 2*(range+(d_state_est./(4*pi/VRBC.lambda)))/3e8;

% Calculate unscaled data vectors.
chirps_per_sym = round(VRBC.Tsym/VRBC.PRI);
samps_per_PRI = round(VRBC.PRI*VRBC.fs);
samps_per_chirp = num_samps;
a = zeros(samps_per_chirp*chirps_per_sym,nnz(M_sym));
for s = 1:nnz(M_sym)
    a(:,s) = (getReturn(tau_tar_est(:,s).', 1, chirps_per_sym, 1/VRBC.fs:1/VRBC.fs:VRBC.Tsym, VRBC.fs, samps_per_PRI, samps_per_chirp, VRBC.beta, VRBC.f0, VRBC.PRI, samps_per_chirp*VRBC.fs));
%     % Bandpass Filter
%     a(:,s) = reshape(bpf_cube(reshape(a(:,s),num_samps,[]), range, VRBC, bpf_width),[],1);
%     % Clutter Filter
%     a(:,s) = reshape(cf(reshape(a(:,s),num_samps,[]), 10, num_chrps),[],1);
end

%% Compare Real/Imag of Modeled to Actual
% Actual Data
seq_delay = VRBC.PRF*2;
data_use = bf_datacube(:,strt_slowtime_ind+seq_delay:strt_slowtime_ind+seq_delay+(VRBC.PRF*VRBC.Tsym*num_sym)-1);
phase_true = phaseEst(data_use, VRBC, range, num_sym*chirps_per_sym, 539, 0, 0);

% Quick View of Micro-Doppler
target.location = range;                         % [m]
window_length = VRBC.Tsym*2;                     % [sec]
window_overlap = 0.8;                            % [ratio]
cube.data = data_use;
cube.config = config;
cube.NFFTrange = num_samps;
cube.NFFTvelocity = num_chrps*num_frames;
cube.ranges = ranges;
cube.dopplers = [-0.5 0.5]*(VRBC.PRF);
STFT.window_len = round(window_length*VRBC.PRF);
STFT.overlap = round(STFT.window_len*window_overlap);
plotDTI(cube, STFT, target)
title(sprintf('Doppler Time Intensity Plot\nReal Data'))

% Modeled Data
prev_sym = zeros(1,N);
state_seq = [1];
prev_sym = zeros(1,N);
a_seq = [a(:,binseq(1)+1)];
d_seq = [d_state_est(:,binseq(1)+1).'];
for s = 2:length(binseq)
    prev_sym = [prev_sym(2:end),binseq(s-1)+1];
    for bit = 1:N
        set_ignore = 1:VRBC.M;
        ign_ind = find(set_N<bit);
        set_ignore= set_ignore(ign_ind);
        if ismember(prev_sym(end-bit+1),set_ignore)==1
            prev_sym(end-bit+1) = 0;
        end
    end
    [~,st_ind] = ismember(prev_sym,states,"rows");
    state_seq = [state_seq, st_ind];
    set_avail = M_sym(st_ind,:);
    ss_ind = nnz(M_sym(1:st_ind-1,:))+nnz(M_sym(st_ind,1:binseq(s)+1));
    a_seq = [a_seq,a(:,ss_ind)];
    d_seq = [d_seq,d_state_est(:,ss_ind).'];
end

% Phase Diff
phase_est = phaseEst(reshape(a_seq, num_samps, []), VRBC, range, num_sym*chirps_per_sym, num_samps, 0, 0);
sym_ph_dif = [];
sym_ph_dif_vec = [];
Lsym = length(a(:,1));
for s = 1:length(binseq)
    ind = (1:chirps_per_sym)+(chirps_per_sym*(s-1));
    sym_ph_dif = [sym_ph_dif, repelem(mean(phase_true(ind))-mean(phase_est(ind)),Lsym)];
    sym_ph_dif_vec = [sym_ph_dif_vec, mean(phase_true(ind))-mean(phase_est(ind))];
end
a_seq = reshape(reshape(a_seq,1,[]).*exp(1j.*sym_ph_dif),Lsym,[]);
phase_est = phaseEst(reshape(a_seq, num_samps, []), VRBC, range, num_sym*chirps_per_sym, num_samps, 0, 0);

% Quick View of Micro-Doppler
target.location = range;                         % [m]
window_length = VRBC.Tsym*2;                     % [sec]
window_overlap = 0.8;                            % [ratio]
cube.data = reshape(a_seq,num_samps,[]);
cube.config = config;
cube.NFFTrange = num_samps;
cube.NFFTvelocity = num_chrps*num_frames;
cube.ranges = ranges;
cube.dopplers = [-0.5 0.5]*(VRBC.PRF);
STFT.window_len = round(window_length*VRBC.PRF);
STFT.overlap = round(STFT.window_len*window_overlap);
plotDTI(cube, STFT, target)
title(sprintf('Doppler Time Intensity Plot\nEstimated Data'))

%% Detection
% Viterbi Trellis Calculations and State Matched Filter Detection
T1 = zeros(S,num_sym);
T2 = zeros(S,num_sym);
mf_state_det_seq = [];
sym_ml_det_seq = [];
for s = 1:num_sym
    y_sym = data_use(:,((s-1)*chirps_per_sym)+1:s*chirps_per_sym);
    y_sym = reshape(y_sym,1,[]).';
    l_stats = zeros(S,1);
    for st = 1:S
        len_inds = find(state_ind_vec==st);
        temp_dif = zeros(length(y_sym),nnz(M_sym(st,:)));
        temp_exp = zeros(nnz(M_sym(st,:)),1);
        temp_inner = 0;
        for m = 1:nnz(M_sym(st,:))
            temp_dif(:,m) = (y_sym-exp(1j.*sym_ph_dif_vec(s)).*a(:,len_inds(m)));
            temp_exp(m) = -(temp_dif(:,m)'*(temp_dif(:,m)));
        end
        for m = 1:nnz(M_sym(st,:))
            temp_inner = temp_inner + exp(temp_exp(m)-max(temp_exp));
        end
        l_stats(st) = max(temp_exp)+log(temp_inner)+log(pis(st,s))-log(nnz(M_sym(st,:)));
    end
    l_stats = real(l_stats);
    l_ind = find(l_stats == max(l_stats));
    mf_state_det_seq = [mf_state_det_seq, l_ind(1)];
    post_prob = zeros(S,1);
    for st = 1:S
        if s > 1
            post_prob(st) = l_stats(st);
            if isnan(post_prob(st))
                post_prob(st) = 0;
            end
            T1(st,s) = max(T1(:,s-1)+log(A(:,st))+post_prob(st));
            tempT2 = find(T1(:,s-1)+log(A(:,st))+post_prob(st)==T1(st,s));
            T2(st,s) = tempT2(1);
        else
            post_prob(st) = l_stats(st);
            if isnan(post_prob(st))
                post_prob(st) = 0;
            end
            T1(st,s) = post_prob(st);
            T2(st,s) = 1;
        end
    end
    temp_exp = zeros(nnz(M_sym(1,:)),1);
    for m = 1:nnz(M_sym(1,:))
        temp = (y_sym-exp(1j.*sym_ph_dif_vec(s)).*a(:,m));
        temp_exp(m) = -(temp'*(temp));
    end
    sym_ml_det_seq(s) = find(temp_exp==max(temp_exp));
end
ind = find(T1(:,end)==max(T1(:,end)));
det_seq = [ind(1)];
for s = 1:num_sym-1
    det_seq = [T2(det_seq(1),end-s+1),det_seq];
end
vit_seq = [];
mfst_seq = [];
for s = 1:num_sym-1
    vit_seq = [vit_seq, states(det_seq(s+1),end)];
    mfst_seq = [mfst_seq, states(mf_state_det_seq(s+1),end)];
end
for s = 1:num_sym-1
    if mfst_seq(s) == 0
        mfst_seq(s) = sym_ml_det_seq(s);
    end
    if vit_seq(s) == 0
        vit_seq(s) = sym_ml_det_seq(s);
    end
end

bin_seq1 = binseq+1;
C1 = confusionmat(bin_seq1(1:end-1),vit_seq)
C2 = confusionmat(bin_seq1(1:end-1),mfst_seq)
C3 = confusionmat(bin_seq1,sym_ml_det_seq)
figure(), imagesc([bin_seq1(1:end-1);vit_seq;mfst_seq;sym_ml_det_seq(1:end-1)])
ax = gca;
ax.YTickLabel = {'','Truth   ','', 'Viterbi ','', 'State ML','','Sym ML',''};
title('Compare Sequence Detection'), colorbar









