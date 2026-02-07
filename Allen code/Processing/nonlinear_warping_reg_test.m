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


%% Calculate the displacements between control point pairs 

cpts.init_disp = cpts.atlas - cpts.data;
% cpts.F = cell(3, 1); % Separate the interpolants in each dimension (here, in 3D)

%% Plot the control points and displacements
figure
scatter3(cpts.data(:, 1), cpts.data(:, 2), cpts.data(:, 3))
hold on
scatter3(cpts.atlas(:, 1), cpts.atlas(:, 2), cpts.atlas(:, 3))

% for ind = 1:cpts.numControlPoints
%     % line([cpts.data(:, 1), cpts.atlas(:, 1)], [cpts.data(:, 2), cpts.atlas(:, 2)], [cpts.data(:, 3), cpts.atlas(:, 3)])
%     line([cpts.data(ind, 1), cpts.atlas(ind, 1)], [cpts.data(ind, 2), cpts.atlas(ind, 2)], [cpts.data(ind, 3), cpts.atlas(ind, 3)])
% end
% hold off

% Plot the original displacement "field"
% figure
quiver3(cpts.data(:, 1), cpts.data(:, 2), cpts.data(:, 3), cpts.init_disp(:, 1), cpts.init_disp(:, 2), cpts.init_disp(:, 3))
hold off
legend('Data points', 'Atlas points', 'Displacement')

%% Create the scatteredInterpolant
% Interpolate between the displacements needed at each data point
% to get a function for the estimated displacement field needed 
% for warping the data to the atlas


% cpts.F = scatteredInterpolant(cpts.data, cpts.init_disp, 'linear');
% cpts.F = scatteredInterpolant(cpts.data, cpts.init_disp, 'natural', 'linear');
% cpts.F = scatteredInterpolant(cpts.data, cpts.init_disp, 'linear', 'none');
cpts.F = scatteredInterpolant(cpts.data, cpts.init_disp, 'linear', 'linear');
% for dir = 1:3
%     cpts.F{dir} = scatteredInterpolant(cpts.data(:, dir), cpts.init_disp(:, dir));
% end

% %%
% test = cpts.F(cpts.data);
% figure;
% quiver3(cpts.data(:, 1), cpts.data(:, 2), cpts.data(:, 3), test(:, 1), test(:, 2), test(:, 3))
%% Apply the scatteredInterpolant to the same grid as the data to get an actual deformation field
ds = size(AA_template_50um); % Data size
[y, x, z] = meshgrid(1:ds(1), 1:ds(2), 1:ds(3));
y = y(:); x = x(:); z = z(:);
tic
D = cpts.F(y(:), x(:), z(:));
toc
testD = reshape(D, ds(1), ds(2), ds(3), 3);
testD(isnan(testD)) = 0;

% %% Plot the post-interpolated displacement field
% figure
% 
% % quiver3(y(:), x(:), z(:), D(:, 1), D(:, 2), D(:, 3))
% testrange = 1:10000;
% quiver3(y(testrange), x(testrange), z(testrange), D(testrange, 1), D(testrange, 2), D(testrange, 3))

%%
testvol = imwarp(PDI_template_rigid_reg, testD, 'cubic');
volumeViewer(testvol)

%% 
test2 = PDI_template_rigid_reg(:);
test3 = reshape(test2, ds(1), ds(2), ds(3));
volumeViewer(test3)