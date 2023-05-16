function [ranges, dopplers, velocities] = computeAxes(axes, NFFTrange, NFFTvelocity)
    ranges = linspace(axes.rangeLims(1), axes.rangeLims(2), NFFTrange);
    dopplers = linspace(axes.dopplerLims(1), axes.dopplerLims(2), NFFTvelocity);
    velocities = linspace(axes.velocityLims(1), axes.velocityLims(2), NFFTvelocity);
end
