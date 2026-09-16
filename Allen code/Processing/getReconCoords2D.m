% Description: output 2D coordinate vectors for reconstructed ultrasound
% data, according to the PData struct that Verasonics uses


function [x_mm, z_mm] = getReconCoords2D(PData, P)
    x_mm = (1:PData.Size(2)) .* PData.PDelta(1) .* P.wl .* 1e3; % x [mm]
    z_mm = (1:PData.Size(1)) .* PData.PDelta(3) .* P.wl .* 1e3; % z [mm]
end