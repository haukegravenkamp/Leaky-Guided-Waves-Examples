function [figC, figA] = plotResults(w,k,att,refCp,refAtt,attThreshold,plotProperties,plotPropertiesRef)

if nargin < 7 
    plotProperties = {};
end
if nargin < 8 
    plotPropertiesRef = {};
end

% default plot properties, can be overwritten by user values
plotPropDefault = {'LineStyle','none','Marker', '.', 'Color', [0.38,0.65,0.76],...
    'MarkerSize',16,'DisplayName','MultiParEig'};
plotPropDefaultRef = {'LineStyle','-','Color','k','LineWidth',1.5};

% attenuation
figA = figure('defaulttextinterpreter','latex','Color','w');
hold all
plot(w/2/pi, att, plotPropDefault{:}, plotProperties{:});
plot(refAtt{1}, refAtt{2}, plotPropDefaultRef{:}, plotPropertiesRef{:});
xlim([0 w(end)/2/pi])
ylim([0 attThreshold])
xlabel('$f$ [MHz]')
ylabel('$\eta$ [dB/m]')
box on

% phase velocities
figC = figure('defaulttextinterpreter','latex','Color','w');
hold all
plot(w/2/pi, w./real(k), plotPropDefault{:}, plotProperties{:});
plot(refCp{1}, refCp{2}, plotPropDefaultRef{:}, plotPropertiesRef{:});
xlim([0 w(end)/2/pi])
ylim([0 10])
xlabel('$f$ [MHz]')
ylabel('$c_p$ [km/s]')
box on

end