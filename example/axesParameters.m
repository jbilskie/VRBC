function axes = axesParameters(config)
    axes.sampleRate = config.profileCfg.outSampleRate*1e3;
    axes.startFreq = config.profileCfg.startFreq*1e9;
    axes.PRI = (config.profileCfg.idleTime+config.profileCfg.rampEndTime)*1e-6;
    axes.PRF = 1/axes.PRI;
    axes.rangeRes = config.c/(2*config.BW*1e6);
    axes.centerFreq = (axes.startFreq)+(((config.profileCfg.rampEndTime*1e-6)*(config.profileCfg.rampSlope*1e12))/2);
    axes.lambdaCenterFreq = config.c/axes.centerFreq;
    axes.rangeLims = [0 (axes.sampleRate*config.c)/(2*config.profileCfg.rampSlope*1e12)];
    axes.dopplerLims = [-0.5 0.5]*axes.PRF;
    axes.velocityLims = (axes.lambdaCenterFreq*axes.PRF/4)*[1 -1];
end

