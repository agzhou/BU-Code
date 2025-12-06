% clear all;
load D:\CODE\Mains\PAT\PATAQParameters_L22_14v.mat;
addpath D:\CODE\DataProcess\BingxueLiu\mPATprocess\SubFunctions;

defaultpath=PATAQInfo.savepath;
[FileName,FilePath]=uigetfile(defaultpath);
load([FilePath, FileName]);

% load extction infor
[Pmeter,delimiterOut]=importdata(['D:\Hb extiction spectra\','Hb spectra_Scott Prahl.txt']);

% load US reference image for generating mask
load([FilePath,'USrefmsk']);
 clear bmsk;
ImgBK = PDIHPdb;
[nx, ny] = size(Xmn(:,:,1));
[nx0, ny0] = size(ImgBK);
[xcoor0, ycoor0] = meshgrid(1:1/ny0:2-1/ny0, 1:1/nx0:2-1/nx0);
[xcoor, ycoor] = meshgrid(1:1/ny:2-1/ny, 1:1/nx:2-1/nx);
ImgBK = interp2(xcoor0, ycoor0, ImgBK, xcoor, ycoor,'same'); % interpolate fluence map;
exist bmsk;
ismsk = ans;
if ismsk==0
    bmsk = selectROI(ImgBK);
end

bmsk = ones(size(ImgBK));

%% select wavelength and import spectra info
prompt={'Selected Wavelength(Example: 700 750 800 850)'};
name='Wavelength for mPAT(nm)';
defaultvalue={num2str([690:10:850])};%{'700 710 760 800 810 860 900'};
numinput=inputdlg(prompt,name,1, defaultvalue);
lam = str2num(numinput{1});

index_lam = find(sum(Pmeter.data(:,1) == lam,2)==1);
eHb = Pmeter.data(index_lam,:);

indmWL = find(ismember(mWL,lam) == 1);
for i = 1: 17
medXmn(:,:,i) = medfilt2(abs(Xmn(:,:,i)));
end
Pm = abs(medXmn(:,:,indmWL)).*bmsk;

%% build extinction matrix: nl * 2n
num_pix = nx*ny;
tic
EHb = [];
EHbR = [];
EHbO = [];
Speye = speye(num_pix);
for i = 1:size(lam,2)
    EHbR = eHb(i,3)*Speye;
    EHbO = eHb(i,2)*Speye;
    EHb = [EHb;EHbR EHbO];
end
toc
Pm = sparse(Pm(:));

%calculate oxygenation level
tic 
cHb = lsqminnorm(EHb,Pm);% cHb1 = EHb\Pm; same: QR solver
toc

% %%% active set, no sparse,calculation time too long
% tic
% [cHb1,resnorm,residual,exitflag,output,lambda] = lsqnonneg(EHb,Pm);
% toc   

%% Apply non-negative constrain
tic
options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter','ConstraintTolerance',1e-9);
lb = sparse(zeros([size(EHb,2),1]));
[cHb,resnorm,residual,exitflag,output,lambda] = lsqlin(EHb,Pm,[],[],[],[],lb,[],[],options);
toc



oxyLevel = calOxyLevel(cHb,nx,ny);
cHbR = reshape(cHb(1:num_pix), nx, ny);
cHbO = reshape(cHb(num_pix+1: end), nx, ny);
cHbT = cHbR + cHbO;

%% figure plot
[VzCmap, VzCmapDn,VzCmapUp, pdiCmapUp, PhtmCmap]=Colormaps_fUS;

figure; 
haxes1 = subplot(121); imagesc(abs(log(squeeze(sum(Xmn,3))/i)).*bmsk,'AlphaData',bmsk); axis image;h=colorbar;%set(get(h,'Title'), 'string','%');
title('Averaged PA distribution');axis off;  colormap(haxes1, parula);%caxis([5 20]);
haxes2 = subplot(122);imagesc(oxyLevel.*bmsk, 'AlphaData',bmsk); axis image;h=colorbar;set(get(h,'Title'), 'string','%');
title('Nonnegative Constrained Oxygenation Map'); axis off; colormap(haxes2, VzCmap);

figure; 
subplot(131); imagesc((log(cHbR)).*bmsk); axis image; colorbar; title('cHbR');
subplot(132); imagesc((log(cHbO)).*bmsk); axis image; colorbar; title('cHbO');
subplot(133); imagesc((log(cHbT)).*bmsk); axis image; colorbar; title('cHbT');

% figure; 
% plot(Pmeter.data(680:900)); hold on; text(Pmeter.data(index_lam),'o','color','r');

%%

figure; 
hAxes1 = axes;
imagesc(ImgBK.*bmsk);
caxis(hAxes1,[(1+0.1*sign(min(ImgBK(:))))*min(ImgBK(:)) max(ImgBK(:))*0.9]);
colormap(hAxes1,gray);
colorbar('TickLabels',{[],[]});
axis equal tight;
axis(hAxes1,'off');

hold on;

% hAxes2 = axes;
% h2 = imagesc(oxyLevel); 
% set(h2,'AlphaData',0.45.*bmsk);%double(bmsk).*abs(oxyLevel-60)>15);%double(bmsk).*abs(oxyLevel)/100);%.*((abs(oxyLevel-100))>10));
% 
% colormap(hAxes2,VzCmap);
% caxis(hAxes2,[0,100])
% colorbar
% hold off
% axis equal tight
% axis(hAxes2,'off')

hAxes2 = axes;
h2 = imagesc(oxyLevel.*(ImgBK>-18).*(oxyLevel>60)); 
set(h2,'AlphaData',0.5.*bmsk.*(ImgBK>-18).*(oxyLevel>60));%double(bmsk).*abs(oxyLevel-60)>15);%double(bmsk).*abs(oxyLevel)/100);%.*((abs(oxyLevel-100))>10));

colormap(hAxes2,VzCmap);
caxis(hAxes2,[60,100])
colorbar
hold off
axis equal tight
axis(hAxes2,'off')

linkaxes([hAxes1,hAxes2]);
%% save data
save([FilePath, 'LinearUnmixingData'],'Pm','cHbR','cHbO','cHbT','oxyLevel');
save([FilePath, 'USrefmsk'],'bmsk','-append');