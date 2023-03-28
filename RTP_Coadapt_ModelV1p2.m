
%--------------------------------------------------------------------------
% This MATLAB code simulates topographic axonal mapping along the anterior/
% posterior axis of the retinotectal system. It is based on the model 
% published in Gebhardt, C., et al., Development, 139, 335 (2012) and 
% modified to include growth cone adaptation. It is associated with the 
% paper Fiedeling, F., et al., „Ephrin-A/EphA specific co-adaptation as a 
% novel mechanism in topographic axon guidance“.

% Copyright (C) 2017, Franco Weth (franco.weth@kit.edu)

% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU general Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version. This program is distributed in the hope that 
% it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
% GNU General Public License for more details.
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Updated Version
% New Features: 
%     - GC and field size adaptable to real life experiments
%     - implementation of a triangular substrate
%     - Quality scaling of substrate independent of fieldsize
%     - option to save GC movement as AVI video
%     - monitoring of different important quantities

clear all

%% USER INPUT--------------------------------------------------------------


% [FILE MANAGEMENT] GENERAL PARAMETERS-------------------------------------

folder          = '04_SmallField_GCBehaviour';
file            = 'GC_5_step_4';
save_Data       = 0;
gif_save        = 0;                                                        % Attention! Can heavily affect performance depending on field size
gif_nthframe    = 10;

S_scaleFactor   = 5;            % for generated substrate resolution      % O(N^2)

                                                                                


% [MEASUREMENT] Scale in micrometer!!
SubstrateSizeX   = 500;
SubstrateSizeY   = 500;
GC_Size          = 5; 
PatternSizeSmall = 2;
PatternSizeLarge = 50;

StepSize = 4;
steps           = 100;         % Iterations of GC growth (ref: 2k)         % O(N)

m_conv_factor = 1/StepSize;

% [SIMULATION] GENERAL PARAMETERS------------------------------------------

SubstrateName = 'Jan_Cones';                                                % 'Substrate_Tectum', 'Substrate_Gap_Assay', 'Substrate_Tectal_Innervation'

FieldSizeX      = ceil(SubstrateSizeX * m_conv_factor);                       % Target field, in which the GC move                                                         
FieldSizeY      = ceil(SubstrateSizeY * m_conv_factor);           

No_GC           = 10;          % Number of GCs used                        % O(N)
sigma_GC        = GC_Size * m_conv_factor;                                 % Standard Deviation of GC distribution



offset          = 1;            % Prevent the GC from getting stuck at edge



% ADAPTATION---------------------------------------------------------------
adap            = 1;            % 1: on
adapHistory     = 10;           % k (number of considered adaptation steps)       
no_adap         = 10;
cis_factor      = 1;
pre_adap        = 1;            
mu              = 0.006;                                                        
lambda          = 0.0045;    
knockIn         = 0;             

% OTHER PARAMETERS (see description)---------------------------------------
GCcutoff        = 0.01;         % Below, Gaussian of GC is assumed 0

ftw = 0;                        % 1: GC forced to take step in x
Qx              = 0;            % [0,1] Driving force                                                      
Qy              = 0; 
sigma_Step      = 0.12;         % Std for Gaussian of step decision
      

x_shift         = 5;            % initial step in x of the GC      
C_dynamic       = 1;            % Boolean, if a dynamic C (below) shall be used
C               = 100;          % Fiber-Fiber factor C(i)                   % Simulates increased terminal number of GC (for FF int.) with iteration

Pedestal_Receptor_Retina    = 0;     
Pedestal_Ligand_Retina      = 0;       
Pedestal_Receptor_Target    = 0;
Pedestal_Ligand_Target      = 0;




o               = 100/FieldSizeX;

kappa_retina    = o*0.025;       % never used?
omega_retina    = 0.4;           % never used?

    

% CREATE NAMES AND FILES---------------------------------------------------
path = fullfile(pwd,'Results', folder);                                    
formatOut = 'yy-mm-dd_HH-MM__';
now_str = datestr(now,formatOut);
file_name = fullfile(path,strcat(now_str,'_',file));

gif_name = strcat(file_name,'_', 'GCR.avi');

% GIF PARAMETERS-----------------------------------------------------------
if gif_save ==1
    factor_plot = 10;
    v = VideoWriter(gif_name, 'Uncompressed AVI');          
    v.FrameRate = 10;    % in fps
    open(v)
end

%% DEFINING VECTORS AND MATRICES-------------------------------------------

FieldSizeXtd    = FieldSizeX + 2*offset;
FieldSizeYtd    = FieldSizeY + 2*offset;

YDrang = [(1-Qy)/3 1/3+(1-Qy)/3 1];                                         % Probabilities for y steps

GrowthConeLigand   = zeros(FieldSizeXtd,FieldSizeYtd);                      % Ligand field of the GC
Lxy   = zeros(FieldSizeXtd,FieldSizeYtd);                                   % Ligand field of the substrate       

GrowthConeReceptor = zeros(FieldSizeXtd,FieldSizeYtd);                      % Receptor field of the GC
Rxy = zeros(FieldSizeXtd,FieldSizeYtd);                                     % Receptor field of the substrate

AxonReceptor = zeros(1,No_GC);                                              % RF
AxonLigand   = zeros(1,No_GC);                                              % LF
AxonReceptor_REF = zeros(1,No_GC);                                          % Reference value 
AxonLigand_REF   = zeros(1,No_GC);

xtHistory  = zeros(steps,No_GC);                                            % Positions of each GC in each iteration
ytHistory  = zeros(steps,No_GC);

DxHistory       = zeros(steps,No_GC);                                       % History values of potential Dx
AbsDxHistory    = zeros(steps,No_GC);                                       % -"- abs value of Dx
QxHistory       = zeros(steps,No_GC);                                       % -"- of forward drive
FactorHistory   = zeros(steps,1);                                           % -"- Factor GC-GC
revHistory = zeros(steps, No_GC);                                           % -"- total rev signal
fwdHistory = zeros(steps, No_GC);                                           % -"- total fwd signal

% ADAPTATION---------------------------------------------------------------
adapmeandenom   = sum(1:adapHistory);

ValAdapRec         = zeros(1,No_GC);
ValAdapLig         = zeros(1,No_GC);
ValAdapRecHistory  = zeros(steps,No_GC);
ValAdapLigHistory  = zeros(steps,No_GC);
AbsAdapRec         = zeros(1,No_GC);
AbsAdapLig         = zeros(1,No_GC);
AbsAdapRecHistory  = zeros(1,No_GC);
AbsAdapLigHistory  = zeros(1,No_GC);

ValResRec          = zeros(steps,No_GC); 
ValResLig          = zeros(steps,No_GC); 
ValResRecHistory   = zeros(steps,No_GC); 
ValResLigHistory   = zeros(steps,No_GC); 

AdapCoeff       = zeros(steps,No_GC);
AdapmeanCoeff   = ones (steps,No_GC);
ResRecCoeff     = ones (steps,No_GC);
ResLigCoeff     = ones (steps,No_GC);

ReceptorHistory = zeros(steps,No_GC);
LigandHistory   = zeros(steps,No_GC);


FSfac= 50/FieldSizeX;                                                       % Factor in the exponent of the retina gradient


%% SUBSTRATES--------------------------------------------------------------
if strcmp('Substrate_Tectum',SubstrateName)==true
    Substrate_Tectum;                                                       % n-t/a-p mapping (Figure 4B)   % original
end
if strcmp('Jan_Cones',SubstrateName)==true
    Jan_Cones;                                                              % 
end    
if strcmp('Substrate_Gap_Assay',SubstrateName)==true
    Substrate_Gap_Assay;                                                    % 
end  
if strcmp('Substrate_Tectal_Innervation',SubstrateName)==true
    Substrate_Tectal_Innervation;                                           % 
end  

%% ALLOCATION OF STARTPOSITIONS--------------------------------------------
% DEFINE POSTITION ON THE RETINA-------------------------------------------
if No_GC == 1                                                               
     YStartPos = ceil(FieldSizeX/2);                                        
else
     YStartPos =  round(linspace(offset,FieldSizeX,No_GC));    
end

% DEFINE GAUSSIAN FUNCTION FOR GCs-----------------------------------------
a=0.5;c=0.5;
[W, V] = meshgrid(1:FieldSizeY, 1:FieldSizeX);   
gaussian=@(x0,y0) exp( - (a*((V-x0)/sigma_GC).^2 ...                   
    + c*((W-y0)/sigma_GC).^2));   

weight=zeros(FieldSizeX, FieldSizeY);
weightxtd=zeros(FieldSizeXtd, FieldSizeYtd);

%% INITIALISATION----------------------------------------------------------
% Define initial positions of the GCs in the target field 
% and the respective potential D
for n=1:No_GC
    % DEFINE INITIAL FIELDS FROM RETINAL POSITION--------------------------
    xt  = 1+x_shift+offset;                                                 % Define Position in Target field xt,yt
                                                                                
    if No_GC == 1                                                    
        yt  = ceil(FieldSizeY/2);
    else
        yt  = round(((FieldSizeY-1)/(No_GC-1))*n  ...
              +((No_GC-FieldSizeY)/(No_GC-1))) + offset;                    % 2. Summand? Freie Felder/ GC
    end
    

    AxonReceptor(n) = (exp(FSfac*0.05*(YStartPos(n)-FieldSizeX/2))+...      % Initialization in Retina according to the exponential
         Pedestal_Receptor_Retina)*pre_adap;
    AxonLigand(n)   = (exp(FSfac*-0.05*(YStartPos(n)-1-FieldSizeX/2))+...
         Pedestal_Receptor_Retina)*pre_adap;
    
%     AxonReceptor(n) = n * 3.4/No_GC;
%     AxonLigand(n) = 3.4 -n * 3.4/No_GC;

    if knockIn > 0                                                          % ???
        for f=1:floor(No_GC/2)
            if n==2*f                                 %Isl EphA3 knockIn
                AxonReceptor(n) = AxonReceptor(n)+knockIn;
                AxonLigand(n)   = 1/AxonReceptor(n);
                break
            else                                      %wildtype
                AxonReceptor(n) = AxonReceptor(n);              
                AxonLigand(n)   = AxonLigand(n);
            end
        end
        
    end


     ReceptorHistory(1,n)= AxonReceptor(n);
     LigandHistory(1,n)= AxonLigand(n);
     
     AxonReceptor_REF(n)= AxonReceptor(n);
     AxonLigand_REF(n)  = AxonLigand(n);
     
   
    % STEP DECISION--------------------------------------------------------
    xrandom = rand;                                                         % Randomly step in xt direction -1,0,+1
    if(xrandom < (1-Qx)/3)                                                  
        xtDirection = -1;
    elseif(xrandom < 1/3+(1-Qx)/3)
        xtDirection = 0;
    else
        xtDirection = 1;
    end
    
    yrandom = rand;                                                         % Randomly step in yt direction -1,0,+1
    if(yrandom < YDrang(1))
        ytDirection = -1;
    elseif(yrandom < YDrang(2))
        ytDirection = 0;
    else
        ytDirection = 1;
    end

    if xt+xtDirection<1+offset                                              % Check if new Position is still in the actual field                 
        xtDirection=0;  
    elseif xt+xtDirection>FieldSizeX+offset
        xtDirection=0;
    end
    if yt+ytDirection<1+offset
        ytDirection=0;
    elseif yt+ytDirection>FieldSizeY+offset
        ytDirection=0;
    end

    Lxy = SubstrateLigand;                                                  % Load Substrate
    Rxy = SubstrateReceptor;
    
    
    Dx1=abs(log(((AxonReceptor(n)*(Lxy(xt,yt)+AxonLigand(n)))/...           % Calculte D for old and new position
        (AxonLigand(n)*(Rxy(xt,yt)+AxonReceptor(n))))));                    % only at specific location xt,yt not surrounding cone
    Dx2=abs(log(((AxonReceptor(n)*(Lxy(xt+xtDirection,yt+ytDirection)+...
        AxonLigand(n)))/(AxonLigand(n)*(Rxy(xt+xtDirection,yt+...
        ytDirection)+AxonReceptor(n))))));
     
    DxHistory(1,n) = Dx1;                                                   
   
    AdapCoeff(1,n) = 1; 
                 
    WDx1 = wkeitpd(Dx1,sigma_Step);                                         % Probabilities for both positions
    WDx2 = wkeitpd(Dx2,sigma_Step);
    
    % MAKE STEP------------------------------------------------------------
    if ftw == 1                                                             % Forced step forwards (see assignment)
        xt=xt+1;        

    elseif ftw == 0                                                         % Random step forwards/sideways (with probabilities)
        if rand>wkeit01(WDx1,WDx2)
            xt = xt+xtDirection;
            yt = yt+ytDirection;
        end
    
    end
    
    xtHistory(1,n) = xt;
    ytHistory(1,n) = yt;
end


%% ITERATION---------------------------------------------------------------
% -------------------------------------------------------------------------

tic
for ii=2:steps

    % Determine C factor   (FF interaction)
    if C_dynamic == 1
        GC_GCfactor=C*(-exp(-log(2.^((ii./(steps./4)).^5)))+1);
    elseif C_dynamic == 0
        GC_GCfactor=C;
    end
        
    GC_SUBfactor=1;
    FactorHistory(ii)=GC_GCfactor;     


    % Calculate total L/R values of ALL Growth cones-----------------------
        % Later used to calculate FF interaction   
    allGrowthConeLigand   = zeros(FieldSizeXtd,FieldSizeYtd);
    allGrowthConeReceptor = zeros(FieldSizeXtd,FieldSizeYtd);

    for nn=1:No_GC                                                   % create allGrowthCone... fields (with individual gcs)
                
        xn=xtHistory(ii-1,nn)-1;                                            % load old target position (minus 1 for gaussian as Matlab sucks and starts indices with 1)
        yn=ytHistory(ii-1,nn)-1;                                            
        weight=gaussian(xn,yn);                                             % Gaussian distribution of concentrations (line 159)                                 
        weight(weight<GCcutoff)=0;
        weightxtd(2:FieldSizeXtd-1,2:FieldSizeYtd-1)= ...
                   weight(1:FieldSizeX,1:FieldSizeY);


        allGrowthConeLigand   = allGrowthConeLigand   + ...
            AxonLigand(nn)  .*weightxtd;
        allGrowthConeReceptor = allGrowthConeReceptor + ...
            AxonReceptor(nn).*weightxtd;
        
    end
    % ---------------------------------------------------------------------
    
    for n=1:No_GC                                                    % Calculate remaining quantities with all GCs
     
        if ii<=no_adap                                                      % fill up R/L history until length of adaptation                                                          
            ReceptorHistory(ii,n)=ReceptorHistory(1,n);                                
            LigandHistory(ii,n)=LigandHistory(1,n);                                        
        end                                                                                

        xt = xtHistory(ii-1,n);
        yt = ytHistory(ii-1,n);       
        
        % DEFINE CONCENTRATIONS OF CURRENT GC------------------------------
        weight=gaussian(xt-1,yt-1);                                         % -1 again because of Matlab
        weight(weight<GCcutoff)=0;
        weightxtd(2:FieldSizeXtd-1,2:FieldSizeYtd-1) = ...
                   weight(1:FieldSizeX,1:FieldSizeY);
        
        currentGrowthConeLigand   = AxonLigand(n)  .*weightxtd;             % Concentration of R/L of current GC (zentraler Wert * Gauss)
        currentGrowthConeReceptor = AxonReceptor(n).*weightxtd;     
        
        meancurrentGCLigand   = sum(sum(currentGrowthConeLigand))/...       % Average Value of current GCL
                                nnz(currentGrowthConeLigand);
        meancurrentGCReceptor = sum(sum(currentGrowthConeReceptor))/...
                                nnz(currentGrowthConeReceptor);     
            
        currentGCLigandRef   = AxonLigand_REF(n)  .*weightxtd;              % Reference Value (1st Axon concentration) for adaptation later         
        currentGCReceptorRef = AxonReceptor_REF(n).*weightxtd;                      %            
                                                                                    %
        meancurrentGCLigandRef   = sum(sum(currentGCLigandRef))/...                 %
                                   nnz(currentGCLigandRef);                         %
        meancurrentGCReceptorRef = sum(sum(currentGCReceptorRef))/...               %
                                   nnz(currentGCReceptorRef);                       %
        
        GrowthConeLigand   = allGrowthConeLigand   - ...
            currentGrowthConeLigand;
        GrowthConeReceptor = allGrowthConeReceptor - ...
            currentGrowthConeReceptor;
        %------------------------------------------------------------------


        % STEP DECISION----------------------------------------------------
        % (same as above)
        xrandom = rand;
        if(xrandom < (1-Qx)/3)
            xtDirection = -1;
        elseif(xrandom < 1/3+(1-Qx)/3)
            xtDirection = 0;
        else
            xtDirection = 1;
        end

        yrandom = rand;
        if(yrandom < YDrang(1))
            ytDirection = -1;
        elseif(yrandom < YDrang(2))
            ytDirection = 0;
        else
            ytDirection = 1;
        end

        if xt+xtDirection<1+offset
            xtDirection=0;
        elseif xt+xtDirection>FieldSizeX+offset
            xtDirection=0;
        end
        if yt+ytDirection<1+offset
            ytDirection=0;
        elseif yt+ytDirection>FieldSizeY+offset
            ytDirection=0;
        end                 
        % -----------------------------------------------------------------

        
        % CALCULATE REV AND FWD VALUES-------------------------------------
        currentGrowthConeLigand_r = imresize( ...                           % Resized images used for increased substrate resolution 
            currentGrowthConeLigand, S_scaleFactor, 'nearest');
        currentGrowthConeReceptor_r = imresize( ...
            currentGrowthConeReceptor, S_scaleFactor, 'nearest');
        GrowthConeReceptor_r = imresize( ...
            GrowthConeReceptor, S_scaleFactor, 'nearest');
        GrowthConeLigand_r = imresize( ...
            GrowthConeLigand, S_scaleFactor, 'nearest');

        rev_FT = GC_SUBfactor.*currentGrowthConeLigand_r.*...               % Calculate individual components of the reverse signal
        SubstrateReceptor;
        rev_FF = GC_GCfactor.*currentGrowthConeLigand_r.*...
        GrowthConeReceptor_r;
        rev_CIS =  cis_factor.*currentGrowthConeLigand_r.*...                      
             currentGrowthConeReceptor_r;

        rev = rev_FT + rev_FF + rev_CIS;                                    % Total reverse signal
%         imshow(rev)
        
        fwd_FT = GC_SUBfactor.*currentGrowthConeReceptor_r.*...             % Calculate individual components of the forward signal
            SubstrateLigand;
        fwd_FF = GC_GCfactor.*currentGrowthConeReceptor_r.*...
            GrowthConeLigand_r;
        fwd_CIS = cis_factor.*currentGrowthConeReceptor_r.*...
             currentGrowthConeLigand_r;   

        fwd = fwd_FT + fwd_FF + fwd_CIS;                                    % Total forward signal
    
        fwdmean=sum(sum(fwd))/nnz(fwd);
        revmean=sum(sum(rev))/nnz(rev);
        
        revHistory(ii,n)=revmean;                                           % Gleich richtig initialisieren!!!
        fwdHistory(ii,n)=fwdmean;                                            
        %------------------------------------------------------------------
        

        % CALCULATE Dx-----------------------------------------------------
        Dx=abs(log(revmean/fwdmean));
          
        DxHistory(ii,n) = Dx;
        % -----------------------------------------------------------------
        

        % Do the same as above for the target field------------------------ % Only used for the adaptation
        weight=gaussian(xt-1+xtDirection,yt-1+ytDirection);
        weight(weight<GCcutoff)=0;
        weightxtd(2:FieldSizeXtd-1,2:FieldSizeYtd-1) = ...
        weight(1:FieldSizeX,1:FieldSizeY);


        currentGrowthConeLigand_target   = AxonLigand(n)  .*weightxtd;
        currentGrowthConeReceptor_target = AxonReceptor(n).*weightxtd;
        
        currentGrowthConeLigand_target_r = imresize( ...
            currentGrowthConeLigand_target, S_scaleFactor, 'nearest');
        currentGrowthConeReceptor_target_r = imresize( ...
            currentGrowthConeReceptor_target, S_scaleFactor, 'nearest');
  

        rev_target = GC_SUBfactor.*currentGrowthConeLigand_target_r.*...
            SubstrateReceptor+GC_GCfactor.*...
            currentGrowthConeLigand_target_r.*GrowthConeReceptor_r+...
            cis_factor.*currentGrowthConeLigand_target_r.*...
            currentGrowthConeReceptor_target_r;
        
        fwd_target = GC_SUBfactor.*currentGrowthConeReceptor_target_r.*...
            SubstrateLigand+GC_GCfactor.*...
            currentGrowthConeReceptor_target_r.*GrowthConeLigand_r+...
             cis_factor.*currentGrowthConeReceptor_target_r.*...
             currentGrowthConeLigand_target_r;
        
        fwdmean_target=sum(sum(fwd_target))/nnz(fwd_target);
        revmean_target=sum(sum(rev_target))/nnz(rev_target);
        
        Dx_target=abs(log(revmean_target/fwdmean_target));

        AbsDxHistory(ii,n) = log(revmean_target/fwdmean_target);
        % -----------------------------------------------------------------
        
%% --------------------------Adaptation------------------------------------

        if adap==1 
            % CALCULATE COEFFICIENTS AND RESETTING FORCE-------------------
            adap_rev = currentGrowthConeLigand_r.*...
                      (currentGrowthConeReceptor_r+SubstrateReceptor);
            adap_fwd = currentGrowthConeReceptor_r.*...
                      (currentGrowthConeLigand_r+SubstrateLigand);
            adap_rev_mean = sum(sum(adap_rev))/nnz(adap_rev);
            adap_fwd_mean = sum(sum(adap_fwd))/nnz(adap_fwd);
            
            AdapCoeff(ii,n)    = 1+log((1+(mu*(abs(log(adap_rev_mean/...
                adap_fwd_mean))))));    
         

            ResRecCoeff(ii,n) = lambda*(meancurrentGCReceptorRef-...
                                meancurrentGCReceptor);     
            ResLigCoeff(ii,n) = lambda*(meancurrentGCLigandRef-...
                                meancurrentGCLigand);

            ResRecCoeff(2,n) = 1;       
            ResLigCoeff(2,n) = 1;
            % -------------------------------------------------------------

                                                                            
            % CALCULATE ADAPTATION COEFF AND UPDATE SENSOR CONCENTRATION    % If Histroy is filled (a and f coefficients) update sensors
            if ii>adapHistory && ii>no_adap                                 
                
                AdapmeanCoeff(ii,n) = 0;
               
                for k=0:adapHistory
                    AdapmeanCoeff(ii,n) = AdapmeanCoeff(ii,n) + ...
                                          k*AdapCoeff(ii-adapHistory+k,n);    
                end
                
               AdapmeanCoeff(ii,n) = AdapmeanCoeff(ii,n)/adapmeandenom;  

               AxonReceptor(n) = AxonReceptor(n)*...
                                 AdapmeanCoeff(ii,n)+ResRecCoeff(ii,n);      
               AxonLigand(n)   = AxonLigand(n)*...
                                 AdapmeanCoeff(ii,n)+ResLigCoeff(ii,n);
               ReceptorHistory(ii,n) = AxonReceptor(n);
               LigandHistory(ii,n)   = AxonLigand(n);
            end
            % -------------------------------------------------------------

% canonical, unproportional adaptation of fwd and rev signals
 
%         AdapCoeff_rec(ii,n)=(log(1/ReceptorHistory(1,n)*...
%                             SubstrateLigand(xt,yt))); 
%         AdapCoeff_lig(ii,n)=(log(1/LigandHistory(1,n)*...
%                             SubstrateReceptor(xt,yt)));
%                         
%         AdapCoeff_rec(ii,n)=(AdapCoeff_rec(ii,n)-(2*AdapCoeff_rec(ii,n)));
%         AdapCoeff_lig(ii,n)=(AdapCoeff_lig(ii,n)-(2*AdapCoeff_lig(ii,n)));
%  
%         
%         AxonReceptor(n)=((AxonReceptor(n)-SubstrateReceptor(xt,yt))*...
%                         exp(-tau*AdapCoeff_rec(ii,n)^2))+...
%                         SubstrateReceptor(xt,yt); 
%         AxonLigand(n)=(AxonLigand(n)-SubstrateLigand(xt,yt))*...
%                        exp(-tau*AdapCoeff_lig(ii,n)^2)+...
%                        SubstrateLigand(xt,yt);   
%                    
%         ReceptorHistory(ii,n) = AxonReceptor(n);
%         LigandHistory(ii,n) = AxonLigand(n);
         
        else
            ReceptorHistory(ii,n) = AxonReceptor(n);
            LigandHistory(ii,n)   = AxonLigand(n);
        
        end

        % MAKE STEP--------------------------------------------------------
        wDx = wkeitpd(Dx,sigma_Step);
        wDx_target = wkeitpd(Dx_target,sigma_Step);

        if ftw == 1
            xt=xt+1;   
            
        elseif ftw == 0
        
            if rand>=wkeit01(wDx,wDx_target)
                xt = xt+xtDirection;
                yt = yt+ytDirection;
            end  
        end

  
        xtHistory(ii,n) = xt;
        ytHistory(ii,n) = yt;
        % -----------------------------------------------------------------
    end

    
   

    %ITERATION FINISHED----------------------------------------------------
    if gif_save == 1                                                        % Update and save the Plot
        if mod(ii,gif_nthframe)==0
            Jan_Update_GCRPlot;
        end   
    end
    
    clc;
    display(steps-ii)

    % PRINT TIMER
    time_elapsed = toc;      
    if ii==gif_nthframe                                                     % Include time to save the avi file in length
          toc1=toc/(gif_nthframe-1); 
    end
    if ii>gif_nthframe                                      
        estimated_time_remaining = (steps-ii)*toc1;
        disp([num2str(round(time_elapsed)),' seconds down, ', ...
            num2str(round(estimated_time_remaining)),' seconds to go!'])
    end
    %----------------------------------------------------------------------
end
if gif_save == 1
    close(v)
end




%% PLOTS

%--------------------------------------------------------------------------
%Map plot


if knockIn > 0      % Isl2-EphA3 constructs
    for i=1:round(No_GC/2)             
        xEndeWildTyp(i,1) = xtHistory(steps, 2*i-1)-offset;
        yStartWildTyp(i,1) = YStartPos(1, 2*i-1)-offset;
        
    end

    for i=1:floor(No_GC/2)
        xEndeKnockIn(i,1) = xtHistory(steps, 2*i)-offset;
        yStartKnockIn(i,1) = YStartPos(1, 2*i)-offset;
        
    end

    hold on
    scatter(xEndeWildTyp, yStartWildTyp, 'b', '*');
    scatter(xEndeKnockIn, yStartKnockIn, 'r', '*');
    axis([0 50 0 50]);
    hold off

elseif knockIn == 0
    
    figure_mapping = figure;
    
    xend=xtHistory(steps,:)-offset;
    xstart=xtHistory(1,:)-offset;
    
    yend=ytHistory(steps,:)-offset;
    ystart=round(YStartPos)-offset;
    
    scatter(xend, ystart, 'b', '*');
    axis([0 FieldSizeX 0 FieldSizeX]);
    title('Mapping of GC Position')
    ylabel('Initial y Position')
    xlabel('Final x Position')
    
    if save_Data == 1
        Jan_SaveData;
    end
end


%--------------------------------------------------------------------------
% % Gap Plot
% blue=[0 0.7 1];
% red =[1 0.3 0.4];
% 
% figure
% hold on
% rectangle('Position',[0,0,unterkante,FieldSizeX],'FaceColor',...
% [1 0.3 0.4]);
% rectangle('Position',[oberkante,0,FieldSizeX-oberkante,FieldSizeX],...
% 'FaceColor',[1 0.3 0.4]);
% scatter(xend,ystart,300,'MarkerEdgeColor','k','MarkerFaceColor',...
% [0.4 0.4 0.4],'LineWidth',1);
% alpha(0.5);
% axis([0 FieldSizeX 0 FieldSizeX]);
% hold off

