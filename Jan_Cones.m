% This Script creates a Receptor and Ligand field with a triangular pattern
% with given large and small ends of a cone

% To generate more accurate results, an ideal triangular pattern is created
% and later resized to the inteded field size by a bicubic interpolation


% USER INPUT---------------------------------------------------------------
cone_size_small = round(PatternSizeSmall * m_conv_factor,1);
cone_size_large = round(PatternSizeLarge * m_conv_factor,1);
%--------------------------------------------------------------------------


cone_period = cone_size_large+cone_size_small;
n_cones = SubstrateSizeY/(PatternSizeSmall+PatternSizeLarge);


% CREATE IDEAL
% PATTERN-----------------------------------------------------% Variables solely related to the pattern use subscript _p and are removed 
p_scaleFactor = 100;                                                       % Create an accurate image with the given parameters, which will later be 
                                                                           % downscaled to the right format
% assign variables and initiate the pattern in pattern
p_cs_small = cone_size_small*p_scaleFactor;                                
p_cs_large = cone_size_large*p_scaleFactor;
p_cone_period = cone_period*p_scaleFactor;

p_x_size = round((p_cs_large-p_cs_small)/2);
p_y_size = round(n_cones * p_cone_period);
pattern = zeros(p_y_size, p_x_size);

y_start=0;

% CREATE THE PATTERN 
for n=0:1:n_cones
    y_start = y_start + p_cs_small/2;                                      % Create black half gap of small cone

    for y1=1:1:(p_cs_large-p_cs_small)/2                                   % Create first triangle of cone
        for xp=1:1:p_x_size
            if xp<y1
                pattern(y1+y_start,xp) = 1;
            end
            if xp==y1
                pattern(y1+y_start,xp) = 0.5;
            end
        end    
    end
    y_start = y_start + (p_cs_large-p_cs_small)/2;                         % Update position


    for y1=1:1:p_cs_small                                                  % Create white gap   
        for xp=1:1:p_x_size
            pattern(y1+y_start,xp) = 1;
        end    
    end
    y_start = y_start + p_cs_small;                                        % Update position
    

    for y1=1:1:(p_cs_large-p_cs_small)/2                                   % Create second triangle of cone
        for xp=1:1:p_x_size
            if xp<y1
                pattern(y_start-y1+(p_cs_large-p_cs_small)/2,xp) = 1;
            end
            if xp==y1
                pattern(y_start-y1+(p_cs_large-p_cs_small)/2,xp) = 0.5;
            end
        end    
    end
    y_start = y_start + (p_cs_large-p_cs_small)/2 + p_cs_small/2;          % Update position
end

pattern = transpose(pattern);
pattern = pattern(1:p_x_size, 1:p_y_size);                                 % Cut out the initial size of the image
% -------------------------------------------------------------------------


% RESHAPE PATTERN TO ACTUAL SIZE-------------------------------------------
SubstrateLigand = imresize(pattern,[S_scaleFactor*FieldSizeXtd, ...
    S_scaleFactor*FieldSizeYtd]);
SubstrateLigand = flipud(SubstrateLigand);
SubstrateReceptor=flipud(SubstrateLigand);                                 % Receptor field is inverse of Ligand
imshow(SubstrateLigand)
% -------------------------------------------------------------------------


% CLEAR UNNECESSARY VARIABLES----------------------------------------------
clear p_scaleFactor p_cs_small p_cs_large p_cone_period 
clear p_x_size p_y_size pattern
% -------------------------------------------------------------------------