% Example matrix
R = transpose(SubstrateReceptor);
L = transpose(SubstrateLigand);

% Convert matrix to RGB image
Mrgb = cat(3, R, zeros(size(R)), L);
Mbw = cat(3, R, R, R);

% Define Plot Parameters
Lines = [1, 9];
% nLines = 3;
% Lines=linspace(1,floor(size(xtHistory,2)-1/nLines),nLines);
nSteps = 100;

% Define the color map based on AxonReceptor_REF
nColors = 256;
cmap = makeColorMap([0 0 1], [1 0 0], nColors);
AxonReceptor_REF = AxonReceptor_REF(:);
colorIndices = ceil((nColors - 1) * AxonReceptor_REF / 3) + 1;

% Set up trace plot figure
tracePlot = figure;
title("GC Trace")
xlabel("xt")
ylabel('yt')
xlim([1, FieldSizeXtd])
ylim([1, FieldSizeYtd])

% Display Mbw as background image
img = Mbw;
image('CData', img, 'XData', [1, FieldSizeXtd], 'YData', [1, FieldSizeYtd])

% Enable hold to overlay plots
hold on;

% Configure colorbar and colormap
ax2 = gca;
clim([0 3]);
colormap(ax2, cmap);
cb2 = colorbar('Location', 'eastoutside');
cb2.Ticks = [0, 1, 2, 3];
cb2.TickLabels = {'0', '1', '2', '3'};
cb2.Label.String = 'Axon Receptor Concentration';
set(ax2, 'ydir', 'normal');

% Define Gaussian function for plotting
[xgc, ygc] = meshgrid(1:FieldSizeYtd, 1:FieldSizeXtd);
gaussian = @(x0, y0) exp(-(a * ((ygc - x0) / sigma_GC) .^ 2 + c * ((xgc - y0) / sigma_GC) .^ 2));

% Plot lines, markers, and Gaussian function
for iLine = Lines
    color = cmap(colorIndices(iLine), :);
    plot(ax2, xtHistory(1:nSteps, iLine), ytHistory(1:nSteps, iLine), '-', 'LineWidth', 2, 'Color', color);
    plot(ax2, xtHistory(1, iLine), ytHistory(1, iLine), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none');
    
    GC_tmp = gaussian(xtHistory(nSteps, iLine), ytHistory(nSteps, iLine));
    GC_tmp(GC_tmp < GCcutoff) = 0;
    GC_tmp_r = imresize(transpose(GC_tmp), 10);
    img_GC = cat(3, color(1) * GC_tmp_r, color(2) * GC_tmp_r, color(3) * GC_tmp_r);
end

% Disable hold and save the figure
hold off;
tracePlot_name = strcat(file_name, '_', 'TracePlot.png');
saveas(tracePlot, tracePlot_name);

% Helper function to make a color map with a smooth gradient between two colors
function cmap = makeColorMap(color1, color2, n)
    cmap = zeros(n, 3);
    for i = 1:3
        cmap(:, i) = linspace(color1(i), color2(i), n);
    end
end
