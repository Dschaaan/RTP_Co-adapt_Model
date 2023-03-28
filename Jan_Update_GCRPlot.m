% % Update the GCR Plot and save to avi movie

img_plot = uint8(10 * imresize(GrowthConeReceptor,[factor_plot*FieldSizeXtd, ...
    factor_plot*FieldSizeYtd]));
set(gca,'nextplot','replacechildren');
plotGCR = imshow(img_plot);
% xlabel('X-Position')
% ylabel('Y-Position')
% title('GCReceptor Concentration')
frame = getframe(gcf);
writeVideo(v,frame);

% Clear unnecessary variables
clear frame

