% Description:
%       Calculate the RX-TX change in plane wave angles for RCA (separate
%       for column-row and row-column)
function delta_theta = calcDeltaThetaRCA(anglesTX, anglesRX)
    naTX = size(anglesTX, 1)/2; % # of transmit angles in one dimension
    naRX = size(anglesRX, 1)/2; % # of receive angles in one dimension
    delta_theta = [];
    % CR
    for tai = 1:naTX
        delta_theta = [delta_theta; anglesRX(1:naRX, :) - anglesTX(tai, :)];
    end
    % RC
    for tai = naTX + 1:naTX*2
        delta_theta = [delta_theta; anglesRX(naRX + 1:naRX*2, :) - anglesTX(tai, :)];
    end
end