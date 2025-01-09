%% function for processing IQ data to BB and IS, ULM
% cluter rejection is based on singular value decomposition (SVD)
function SCC_BBISC(datapath, filename)
% IQ: IQ data
% PRMT: data processing parameter
% vUS: Obtained vUS results
load([datapath,'BBISC-PRSinfo.mat'])
load([datapath,filename]);
%%
[nz,nx,nt]=size(BB);
zShiftInAll=zShift+repmat(zBLKshift(1,:,:),[nt,1,1]);
xShiftInAll=xShift+repmat(xBLKshift(1,:,:),[nt,1,1]);
tNow=datetime('now');
disp(['Processing: ',filename,', ', datestr(now, 'HH:MM:SS')])
ImgRfn=PRSinfo.BBlPix/PRSinfo.BBIScPix;
% BBIS=zeros(nz*ImgRfn,nx*ImgRfn);
BBISC=uint8(zeros(nz*ImgRfn,nx*ImgRfn));
for it=1:nt
    BBIS0=zeros(nz*ImgRfn,nx*ImgRfn);
    itBB=imresize(squeeze(BB(:,:,it)),ImgRfn);
    izShift=round(imresize(squeeze(zShiftInAll(it,:,:)),[nz*ImgRfn,nx*ImgRfn])/PRSinfo.BBIScPix);
    ixShift=round(imresize(squeeze(xShiftInAll(it,:,:)),[nz*ImgRfn,nx*ImgRfn])/PRSinfo.BBIScPix);
    [zB,xB]=find(itBB>0);
    nBB=size(zB,1);
    for iBB=1:nBB
        zB(iBB)=zB(iBB)-izShift(zB(iBB),xB(iBB));
        xB(iBB)=xB(iBB)-ixShift(zB(iBB),xB(iBB));
        BBIS0(zB(iBB),xB(iBB))=1;
    end
    BBISC(:,:)=BBISC+uint8(BBIS0);
end

%%
pathInfo=strsplit(datapath,'/');
fileInfo=strsplit(filename,'-');
SavePath=['/',strjoin(pathInfo(1:end-2),'/'),'/B-BBISC/'];
if ~exist(SavePath)
    mkdir(SavePath);
end
P=PRSinfo;
SaveName=['BBISC-',strjoin(fileInfo(2:end),'-')];
save([SavePath,SaveName],'-v7.3','BBISC','BB','zShift','xShift','zBLKshift','xBLKshift','P');
disp([SaveName]);

