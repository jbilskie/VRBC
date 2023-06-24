function out = simulateCltReturn(radar_vel, radar_vib, var_clts, range_clts, az_angs_clts, full_ts, beta, f0, PRI, T, array_pos)
% Author: Jessica Centers - June 20, 2023
% Description: Simulate FMCW clutter return given FMCW waveform parameters;
%   power, ranges, and angles of multiple clutter points; a description of 
%   any bulk radar movement and vibrations, and horizontal linear recieve 
%   array element positions.
% Inputs:
%   radar_vel: 2 x 1 vector the velocity magntiude and then direction of the radar platform [meters/sec, degrees azimuth]
%   radar_vib: time series of the radar platform vibrations, assuming "forward" direction [meters]
%   var_clts: vector of the power of each clutter point being modeled
%   range_clts: vector of the beginning ranges of each clutter point being modeled [meters]
%   az_angs_clts: vector of the azimuthal angles of each clutter point being modeled [degrees] 
%   full_ts: full fast-time samples vector [sec]
%   beta: chirps slope [Hz/sec]
%   f0: starting frequency of the chirp [Hz]
%   PRI: pulse repetition interval [sec]
%   T: active duration of the PRI [sec]
%   array_pos: vector of element positions [+ is right of the origin] of the array [m]
% Output:
%   out: radar cube (element x fast-time x slow-time) of simulated clutter data

    BW = beta*T;                    % chirp bandwidth [Hz]
    fs = round(1./(full_ts(2:end)-full_ts(1:end-1))); 
    if length(unique(fs))~=1
        warning('The slow-time time series provided is not uniformly sampled.')
    end
    fs = mean(fs);                  % ADC sampling freuency [Hz]
    num_chirps = (max(full_ts)-min(full_ts)+(1/fs))/PRI;
    samps_per_PRI = round(PRI*fs);  % samples per PRI [#]
    if mod((PRI*fs),1)~=0
        warning('PRI is not perfectly divisible by the sampling frequency indicated by the slow-time time series provided.')
    end
    samps_per_chirp = round(T*fs); % samples per chirp [#]
    if mod((T*fs),1)~=0
        warning('Active chirp duration is not perfectly divisible by the sampling frequency indicated by the slow-time time series provided.')
    end
    if mean(array_pos)~=0
        warning('Array does not appear to be centered at the assumed origin.')
    end
    num_clt_pts = length(var_clts);
    if length(unique([length(var_clts),length(range_clts),length(az_angs_clts)]))~=1
        warning('Clutter point powers, ranges, and angles do not indicate the same number of clutter points.')
    end
    
    out = zeros(length(array_pos),samps_per_chirp,num_chirps);
    for cpt = 1:num_clt_pts
        clt_phse = 2*pi*rand(1);
        a = sqrt(var_clts(cpt))*exp(1j*clt_phse);
        for elem = 1:length(array_pos)
            % Provide delay due to element position.
            r = range_clts(cpt)+0.5*(array_pos(elem)*sind(az_angs_clts(cpt)));
            % Calculate the range over time given velocity of the radar and radar platform vibration
            r = r-((radar_vel(1).*full_ts).*cosd(az_angs_clts(cpt)+radar_vel(2)))-radar_vib.*cosd(az_angs_clts(cpt));
            % Delay
            tau = 2.*r./3e8;
            % Calculate waveform
            for g = 0:num_chirps-1
                tau_temp = tau(samps_per_PRI*g+(1:samps_per_chirp));
                t_temp = full_ts(samps_per_PRI*g+(1:samps_per_chirp));
                out(elem,:,g+1) = out(elem,:,g+1)+...
                                  a.*exp(1j*2*pi*((f0.*tau_temp)+(beta.*tau_temp.*(t_temp-(g*PRI)))-(0.5*beta.*tau_temp.^2)));
            end
        end
    end
end