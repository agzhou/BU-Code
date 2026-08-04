%% correlation coefficient calculation, CPU
% DAT: 2D matrix, [nVox,nt]
% DAT2: [1,nt]
function R= Coef(DAT1, DAT2)
[nVox,nt]=size(DAT1);
Numerator=nt.*(sum(DAT1.*DAT2,2))-sum(DAT1,2).*sum(DAT2,2);
Denomenator=sqrt((nt*sum(DAT1.^2,2)-(sum(DAT1,2)).^2).*(nt*sum(DAT2.^2,2)-(sum(DAT2,2)).^2));
R=Numerator./Denomenator;
