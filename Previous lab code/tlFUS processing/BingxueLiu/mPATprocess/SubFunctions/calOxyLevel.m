function oxyLevel = calOxyLevel(cHb,nz,nx)
num_pix = size(cHb,1)/2;
cHbR = cHb(1: num_pix); cHbO = cHb(num_pix+1: end);
oxyLevel = 1./(1+cHbR./cHbO)*100;
oxyLevel = reshape(oxyLevel, [nz,nx]);
