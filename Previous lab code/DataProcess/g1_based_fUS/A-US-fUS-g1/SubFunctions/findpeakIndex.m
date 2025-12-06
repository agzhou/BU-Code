function [locmin,Index,dir] = findpeakIndex(iGG, peakrange)
% max imag for pi/2 
% min imag for pi3/2 
% min real for pi
switch numel(size(iGG))
    case 2 % single pixel
dir = sum(mean(iGG(1:5,:),2));
if dir>0
    dir = 1;
else
    dir= -1;
end
iGG_bar = iGG.*dir;
tic
[locmin,Index] = max(iGG_bar(peakrange(1):peakrange(2),:),[],1);
toc
    case 4 % whole image
tic
dir = sum(mean(iGG(:,:,1:10,:),4),3);
dir(dir>0) = 1;
dir(dir<0) = -1;
iGG_bar = iGG.*dir;
[locmin,Index] = max(iGG_bar(:,:,peakrange(1):peakrange(2),:),[],3);
toc
end

% igg = iGG;
% dir = sum(mean(iGG(1:5,:),2));
% nRpt = size(iGG,2);
% %tic
% if dir > 0
%  for k = 1: nRpt
%    [pks,locs,w,p] = findpeaks(-igg(:,k));
% %    Index(:,k) = locs(find(p == max(p)));
%    Index(:,k) = locs(1);
%     locmin(:,k) = pks(1);
%  end
% else
%  for k = 1: nRpt
%    [pks,locs,w,p] = findpeaks(igg(:,k));
% %    Index(:,k) = locs(find(p == max(p)));
%    Index(:,k) = locs(1);
%    locmin(:,k) = pks(1);
%  end
% end
% %toc