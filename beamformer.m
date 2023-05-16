function W = beamformer(N, axes, data, theta, f)
    lambda = (3e8/f);
    % Calculate distance between antenna elements
    if N==1
        d=0;
    elseif N==2
        d=[-lambda/4, lambda/4];
    elseif N==3
        d=[-lambda/2,0,lambda/2];
    elseif N==4
        d=(-3*lambda/4:lambda/2:3*lambda/4);
    else
        display('Number of receive elements must be between 1 and 4.')
        keyboard
    end
    dd = d(2)-d(1);
    K0 = 2*pi*dd*sind(theta)/lambda;
    angs = -180:1:180;
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