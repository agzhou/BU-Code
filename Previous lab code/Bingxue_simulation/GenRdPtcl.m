function [xp,zp] = GenRdPtcl(xl, zl, pixsize, density, agl)
%%%% last modified by Bingxue liu
agl = agl/180*pi;

xCoor = (-xl/2:pixsize:xl/2);
zCoor = (-zl/2:pixsize:zl/2);
[Xp, Zp] = meshgrid(xCoor, zCoor);

nx = length(xCoor);
nz = length(zCoor);

np = round(xl*zl*density);

temp = rand(nx*nz,1);
dtemp = sort(temp,'descend');
ind = find(temp>dtemp(np+1));

xp0 = Xp(ind);
zp0 = Zp(ind);
rotm = [cos(agl), -sin(agl); sin(agl), cos(agl)];
ptcl = rotm*[xp0';zp0'];
xp = ptcl(1,:);
zp = ptcl(2,:);
end

