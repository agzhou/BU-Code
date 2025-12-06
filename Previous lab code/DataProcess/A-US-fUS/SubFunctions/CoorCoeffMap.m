%% correlation coefficient and map plot
function actMap=CoorCoeffMap(data, Stim, CoefPlot)
if nargin<3
    CoefPlot=0;
end
StimN=squeeze(Stim-mean(Stim));
StimN=StimN./sqrt(sum(StimN.^2));

% cData0=abs(data);
% cData0=movmean(abs(data),3,3);
% cData0(cData0<0.5)=0;
cData0=data;
% cData0(cData0<mean(cData0(:))-std(cData0(:)))=0;
% [nz,nx,nt]=size(cData0);
% cData=movmean(cData0,3,3);
cData=cData0;
cDataN=cData-mean(cData,3);
cDataN=cDataN./sqrt(sum(cDataN.^2,3));

nt = length(StimN);
Start=1;
End=length(StimN);
coefMap=sum(cDataN(:,:,Start:End).*permute(StimN(:,Start:End),[1 3 2]),3);
coefMap(find(isnan(coefMap)==1))=0;
medcoefMap = medfilter(coefMap,5);
z=sqrt(nt-3)*log((1+medcoefMap)./(1-medcoefMap))/2;
zMsk=zeros(size(z));
zMsk(z>3.1)=1; % 1.6

% coefMap_af = medfilter(coefMap,5);

%%%% modified by Bingxue Liu
% medcoefMap = ordfilt2(coefMap,5,ones(3,3)); % median filter 

conFzMsk = bwareafilt(logical(zMsk),[20,inf]);%10
conFzMsk0 = logical(zMsk)-conFzMsk;  % 9 pixels connectivity filter 
actMap=medcoefMap.*conFzMsk;
BB=ones(7,7);
BB(4,4)=49;
BB=BB/49;
actMap_covn = convn(actMap,BB,'same');

figure; imagesc(z);axis image;colorbar;colormap(jet); caxis([-16,16]); title('Z-Score Map');
figure; imagesc(zMsk);axis image; colorbar; title('Z-Score>3.1');
figure; imagesc(medcoefMap);axis image; colorbar;colormap(jet);caxis([-1,1]);title('Correlation Map');
figure; imagesc(actMap);axis image; colorbar; title('Masked Correlation Map(Z-Score>3.1)');caxis([0,1]); 
% 
% %%%% end
% actMap = actMap_covn;
if CoefPlot==1
    %% figure plot
    [VzCmap, VzCmapDn,VzCmapUp]=Colormaps_fUS;
    figure;
    hAxes1 = axes;
    if max(abs(data(:)))>1e6
        ImgShow=abs(mean(data,3)).^0.25;
    else
        ImgShow=abs(mean(data,3));
    end
    imagesc(ImgShow);
    caxis(hAxes1,[(1+0.2*sign(min(ImgShow(:))))*min(ImgShow(:)) max(ImgShow(:))*0.7]);
    colormap(hAxes1,gray)
    colorbar('Ticks',[-10 10], 'TickLabels',{[],[]} );
    axis equal tight
    
    BKmsk=zeros(size(ImgShow));
    BKmsk(abs(ImgShow)>2.5)=1;
    hold on;
    hAxes2 = axes;
    h2=imagesc(actMap);
    set(h2,'AlphaData',abs(actMap)/(max(abs(actMap(:)))/4).*((abs(actMap))>0.25).*BKmsk)
    colormap(hAxes2,VzCmap);
    caxis(hAxes2,[-0.6 0.6])
    % colormap(hot);
    % caxis(hAxes2,[0.2 0.6])
    colorbar
    hold off
    axis equal tight
    axis(hAxes2,'off')
    linkaxes([hAxes1,hAxes2])
end
end

function af=medfilter(a,n)
[nz,nx]=size(a);
af=a;
for iz=1+n:nz-n
    for ix=1+n:nx-n
        tmp=a([-n:n]+iz,[-n:n]+ix);
        st=std(tmp(:));
        m=median(tmp(:));   
        if abs(af(iz,ix)-m)>1.5*st
        af(iz,ix)=m;
        end
    end
end
end
