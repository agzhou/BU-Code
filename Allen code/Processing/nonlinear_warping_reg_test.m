% Nonlinear warping registration

%% Choose control point pairs from the rigidly-registered data and the atlas
US_Atlas_Reg % (Use the orthosliceviewer button)

%% Extract control point pairs
cpts.numControlPoints = size(controlPoints, 1);
% Initialize matrices for the control points. Dimensions [# control point pairs, 3 (xyz)]
cpts.data = NaN(cpts.numControlPoints, 3); % Coordinates from the data to be warped
cpts.atlas = NaN(cpts.numControlPoints, 3); % Corresponding coordinates from the atlas
dpcn = 2; % Data points column #
apcn = 3; % Atlas points column #

for ind = 1:cpts.numControlPoints
    % cpts.data = controlPoints(:, dpcn);  % Coordinates from the data to be warped
    % cpts.atlas = controlPoints(:, apcn); % Corresponding coordinates from the atlas
    
    datapts_split_temp = split(controlPoints{ind, dpcn}, [',', ' ']);
    atlaspts_split_temp = split(controlPoints{ind, apcn}, [',', ' ']);

    for dir = 1:3 % Separate the xyz since they are each a cell from the split
        cpts.data(ind, dir) = str2double(datapts_split_temp{dir});
        cpts.atlas(ind, dir) = str2double(atlaspts_split_temp{dir});
    end
end
clearvars ind dir datapts_split_temp atlaspts_split_temp

%% Calculate the displacements between control point pairs and create the scatteredInterpolant
cpts.init_disp = cpts.atlas - cpts.data;
% cpts.F = cell(3, 1); % Separate the interpolants in each dimension. 3D

% Interpolate between the displacements needed at each data point
% to get a function for the estimated displacement field needed 
% for warping the data to the atlas


cpts.F = scatteredInterpolant(cpts.data, cpts.init_disp, 'linear');
% for dir = 1:3
%     cpts.F{dir} = scatteredInterpolant(cpts.data(:, dir), cpts.init_disp(:, dir));
% end

%% Apply the scatteredInterpolant to the same grid as the data to get an actual deformation field
ds = size(AA_template_50um); % Data size
[y, x, z] = meshgrid(1:ds(1), 1:ds(2), 1:ds(3));
D = cpts.F(y, x, z);
testD = reshape(D, ds(1), ds(2), ds(3), 3);

%%
testvol = imwarp(PDI_template_rigid_reg, testD);