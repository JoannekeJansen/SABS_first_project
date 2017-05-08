close all
clear all
nLines = 18;
legend_str = cell(nLines,1);

myColors = jet(nLines);
t = 0:pi/64:3*pi;
dPhi = pi/16;
lineStyles = {'-', '--', ':', '-.'};
for ii=1:nLines,
    plot(t,sin(t+dPhi*ii),...
        'linestyle', lineStyles{rem(ii-1,numel(lineStyles))+1},...
        'color', myColors(ii,:),...
        'linewidth',3);
    hold on;
    legend_str{ii} = num2str(ii);
end
axis([0 3*pi -1.15 1.6])
legend(legend_str,'location','NorthWest')
for ii=1:nLines,
    plot(t,sin(t+dPhi*ii),...
        'linestyle', lineStyles{rem(ii-1,numel(lineStyles))+1},...
        'color', myColors(ii,:),...
        'linewidth',3);
    hold on;
    legend_str{ii} = num2str(ii);
end
axis([0 3*pi -1.15 1.6])
columnlegend(6,legend_str,'NorthWest');
shg
