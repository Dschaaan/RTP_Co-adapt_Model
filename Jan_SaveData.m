% Save all necessary files and plots

plot_name = strcat(file_name,'_', 'Mapping.png');
gv_name = strcat(file_name, '_', 'General_Variables');
sub_name = strcat(file_name, '_', SubstrateName,'.png');




save(gv_name, 'No_GC','GC_Size','steps','FieldSizeX',...
    'FieldSizeY','GCcutoff','Qy','sigma_Step','mu','lambda','knockIn',...
    'cis_factor','pre_adap','no_adap','x_shift','C_dynamic','C',...
    'Pedestal_Receptor_Retina','Pedestal_Ligand_Retina',...
    'Pedestal_Receptor_Target','Pedestal_Ligand_Target',...
    'omega_retina','ftw','adap','adapHistory', 'SubstrateName')
save(file_name)
saveas(figure_mapping, plot_name);

figure_substrate = figure;
imshow(SubstrateLigand)
title('Substrate Ligand')
saveas(figure_substrate, sub_name)