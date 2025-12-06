function Imgbar = jitCor(Img)
%% post-processing: remove jitter, for laser master PAT data 
[nr, nc, nt] = size(Img);
Imgk = zeros(nr+1,nc,nt);
Imgp = Img;
[~,MaxP] = find(Img(:,:,1)==max(max(Img(:,:,1)))); % column position of maximum value of the first image
for i = 1: nt
    indMax(i) = find(Img(:,MaxP,i)==max(Img(:,MaxP,i)));
    MaxMost = mode(indMax);
    if indMax(i)> MaxMost
%        k = indMax(i)-MaxMost;
%        Imgk(:,:,i) = padarray(Img(:,:,i), k, 0, 'post');
        Imgk(:,:,i) = vertcat(Img(:,:,i),Img(nr, :, 1));
        Imgp(:,:,i) = Imgk(2:end,:,i);
    elseif indMax(i) < MaxMost
        Imgk(:,:,i) = vertcat(Img(1, :, 1),Img(:,:,i));
        Imgp(:,:,i) = Imgk(1:end-1,:,i);
    end
end
Imgbar = sum(Imgp, 3)/nt;
Imgvar = std(Imgp, 0, 3);
Imgsnr = 20*log10(Imgbar./Imgvar);
%% plots
% for i = 1: nt
%     figure(1);
%     plot(Img(:,MaxP,i));
%     hold on;
%     figure(2);
%     plot(Imgp(:,MaxP,i));
%     hold on;
% end
% figure; imagesc(Imgbar.^0.4); axis image; colorbar; title('corrected averaged image');
% figure; imagesc(mean(Img,3).^0.4); axis image; colorbar;title('original averaged image');
% figure; imagesc(Imgvar); axis image; colorbar; title('variance image');
% figure; imagesc(Imgsnr); axis image; colorbar; title('snr image');
end