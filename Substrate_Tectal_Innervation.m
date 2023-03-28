%--------------------------------------------------------------------------
% This MATLAB code contains the target matrix called by 
% RTP_Coadapt_ModelV1p2.m simulating the a-p graded distributions of EphAs 
% and ephrin-As on the tectum with an additional, cue-free stretch in front 
% of the tectums’ anterior pole.

% Copyright (C) 2017, Franco Weth (franco.weth@kit.edu)

% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU general Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version. This program is distributed in the hope that 
% it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
% GNU General Public License for more details.
%--------------------------------------------------------------------------

%% Target Substrate with free space in front of RT/LT counter-gradients

SubstrateLigand=zeros(FieldSizeXtd,FieldSizeYtd);
SubstrateReceptor=zeros(FieldSizeXtd,FieldSizeYtd);

b=o*0.03;


for zi=21:1:FieldSizeX
    zz=zi-19;
    for zn=1:1:FieldSizeY
        SubstrateLigand(zi+offset,zn+offset) = exp(b*(zi-19-102/2));
        SubstrateReceptor(zi+offset,zn+offset) = exp(-b*(zi-20-102/2));
    end
end
imshow(SubstrateReceptor)
