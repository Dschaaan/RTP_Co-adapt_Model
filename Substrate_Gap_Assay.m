%--------------------------------------------------------------------------
% This MATLAB code contains the target matrix called by 
% RTP_Coadapt_ModelV1p2 simulating various ephrin-A/EphA single and 
% double-cue in-vitro gap substrates.

% Copyright (C) 2017, Franco Weth (franco.weth@kit.edu)

% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU general Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version. This program is distributed in the hope that 
% it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
% GNU General Public License for more details.
%--------------------------------------------------------------------------

%% Target Substrate with homogeneous distributions of RT/LT with a free gap

SubstrateLigand=zeros(FieldSizeXtd,FieldSizeYtd);
SubstrateReceptor=zeros(FieldSizeXtd,FieldSizeYtd);

r=0;    %RT
l=4;    %LT


lueckenbreite=40;       
unterkante=ceil(0.4*FieldSizeX);            
oberkante=unterkante+lueckenbreite;


for zi=1:1:FieldSizeY
    for zn=1:1:unterkante
        SubstrateLigand(zn+offset,zi+offset) = l;
        SubstrateReceptor(zn+offset,zi+offset) = r;
    end
        
    for zn=oberkante+1:1:FieldSizeX
        SubstrateLigand(zn+offset,zi+offset) = l;
        SubstrateReceptor(zn+offset,zi+offset) = r;
    end
end
imshow(SubstrateLigand)
