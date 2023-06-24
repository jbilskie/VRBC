function W = beamformer(d, theta, lambda)
    dd = d(2)-d(1);
    K0 = 2*pi*dd*sind(theta)/lambda;
    angs = -180:1:180;
    N = length(d);
    W = ones(N,length(angs));
    for ii = 1:length(angs)
        phase_diff = 2*pi*dd*sind(angs(ii))/lambda;
        for i = 1:N
            W(i,ii) = W(i,ii)*conj(exp(1j*(i-1)*K0))*exp(1j*(i-1)*phase_diff);
        end
    end
    A = sum(W,1);
    G = 10*log((abs(A).^2)/(N^2));
    ind = find(angs==0);
    W = W(:,ind);
end