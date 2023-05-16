function out = getReturn(tau, alpha, num_chirps, time, fs, samps_per_PRI, samps_per_chirp, beta, f0, PRI, T)
    out = [];
    for g = 0:num_chirps-1
        delay_samps = ceil(fs*tau((g*samps_per_PRI)+1));
        chirp_tau = tau((g*samps_per_PRI)+1+delay_samps:(g*samps_per_PRI)+samps_per_chirp);
        chirp_time = time(1+delay_samps:samps_per_chirp);
        if delay_samps > 0
            out = [out, zeros(1,delay_samps), ...
                        exp(-1j*2*pi*((0.5*beta.*(chirp_tau.^2))-(beta.*(chirp_time)+f0).*chirp_tau)), ...
                        zeros(1,round(fs*(PRI-T)))];
        else
            out = [out,  ...
                         exp(-1j*2*pi*((0.5*beta.*(chirp_tau.^2))-(beta.*(chirp_time)+f0).*chirp_tau)), ...
                         zeros(1,round(fs*(PRI-T)))];
        end
    end
    out = alpha.*out;
end