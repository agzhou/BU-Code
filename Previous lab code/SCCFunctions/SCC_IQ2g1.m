% calculate 3D g1, return nz_nx_ny_ntau
% nt: number of time points for ACF calculation, usuall, nt=nxRpt
function SCC_IQ2g1(datapath, filename, g1Info)
load([datapath,filename]);
[nz, nx, nAgle,nRpt] = size(IQ) ;
cIQ=hilbert(squeeze(sum(IQ,3))); % Coherence Compounding IQ data, then Hilbert transform to get complext Coherence compounded IQ data (envelop+phase)
tStart=g1Info(1);
nTau=g1Info(2);
nt=g1Info(3);
% Collect ACFs at each lag i
%%%%%%%%%%%%%%%%%%%
%%%% Jonghwan ACF mean((yi*)*(y(i+tau))/mean(yi**yi);
temp=zeros(nz,nx,nt);
Numer=zeros(nz,nx,nTau);
fileinfo=strsplit(filename,'-');
for itt=tStart:tStart-1+nt
    temp(:,:,itt-tStart+1)=squeeze(cIQ(:,:,itt))'.'.*squeeze(cIQ(:,:,itt));
end
Denom=mean(temp,3); % calculate the denomenator
for it = 1:nTau
    %% calcualte numerator 
    temp=zeros(nz,nx,nt);
        for itt=tStart:tStart-1+nt-it
            
            temp(:,:,itt-tStart+1)=squeeze(cIQ(:,:,itt)'.').*squeeze(cIQ(:,:,it+itt));
        end
    Numer(:,:,it)=sum(temp,3)/(nt-it);
    GG(:,:,it) = Numer(:,:,it)./Denom; 
end
GGsaveName=[strjoin(fileinfo(1:10),'-'),'-',num2str(nTau),'-',num2str(nt),'-GG.mat'];
save([datapath,GGsaveName],'GG','P','NA','g1Info','xCoor','zCoor');
disp('GG saved!')
% %%%% Traditional abs ACF sum(abs(yi-y_mean)*abs(y(i+tau)-y_mean))/sum(yi-t).^2;
% Denom=sum(abs((RR_iy(:,:,:)-repmat(mean_RR,[1 1 nt]))).^2,3); % calculate the denomenator
% for it = 1:ntau
%     %% calcualte numerator 
%     Numer(:,:,it)=sum(abs(RR_iy(:,:,1:nt-it)-repmat(mean_RR,[1 1 nt-it])).*abs(RR_iy(:,:,it+1:nt)-repmat(mean_RR,[1 1 nt-it])),3);
%     GG(:,:,it) = Numer(:,:,it)./Denom; 
% end
%%%%%%%%%%%%%%%%%%%%%%
% %%%% Traditional ACF sum((yi-y_mean)*(y(i+tau)-y_mean))/sum(yi-t).^2;
% Denom=mean(abs((RR_iy(:,:,:)-repmat(mean_RR,[1 1 nt]))).^2,3); % calculate the denomenator
% for it = 1:ntau
%     %% calcualte numerator 
%     Numer(:,:,it)=sum((RR_iy(:,:,1:nt-it)-repmat(mean_RR,[1 1 nt-it])).*(RR_iy(:,:,it+1:nt)-repmat(mean_RR,[1 1 nt-it])),3)/(nt-it);
%     GG(:,:,it) = Numer(:,:,it)./Denom; 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% modify Jonghwan ACF (yi-mean(yi));
% temp=zeros(nz,nx,nt);
% Numer=zeros(nz,nx,ntau);
% for itt=1:nt
%     temp(:,:,itt)=(squeeze(RR_iy(:,:,itt))'.'-mean_RR).*(squeeze(RR_iy(:,:,itt))-mean_RR);
% end
% Denom=mean(temp,3); % calculate the denomenator
% for it = 1:ntau
%     %% calcualte numerator 
%     temp=zeros(nz,nx,nt);
%     for itt=1:nt-it
%         temp(:,:,itt)=(squeeze(RR_iy(:,:,itt)'.')-mean_RR).*(squeeze(RR_iy(:,:,it+itt))-mean_RR);
%     end
%     Numer(:,:,it)=sum(temp,3)/(nt-it);
%     GG(:,:,it) = Numer(:,:,it)./Denom; 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Tradition ACF sum((yi-y_mean)*(y(i+tau)-y_mean))/sum(yi-t).^2;
% Denom=sum((abs(RR_iy(:,:,:)-repmat(mean_RR,[1 1 nt]))).^2,3); % calculate the denomenator
% for it = 1:ntau
%     %% calcualte numerator 
%     Numer(:,:,it)=sum((RR_iy(:,:,1:nt-it)-repmat(mean_RR,[1 1 nt-it])).*(RR_iy(:,:,it+1:nt)-repmat(mean_RR,[1 1 nt-it])),3);
%     GG(:,:,it) = Numer(:,:,it)./Denom; 
% end
%%%%%%%%%%%
% for it = 1:ntau
%     %% calcualte numerator 
%     Numer(:,:,it)=sum((RR_iy(:,:,1:nt-it)-repmat(mean_RR,[1 1 nt-it])).*(RR_iy(:,:,it+1:nt)-repmat(mean_RR,[1 1 nt-it])),3);
% %     Numer(:,:,it)=sum((RR_iy(:,:,1:nt-it)-repmat(mean_RR,[1 1 nt-it])).*(RR_iy(:,:,it+1:nt)-repmat(mean_RR,[1 1 nt-it])),3)/(nt-it);
% %     Numer(:,:,it)=mean((RR_iy(:,:,1:nt-it)-repmat(mean_RR,[1 1 nt-it])).*(RR_iy(:,:,it+1:nt)-repmat(mean_RR,[1 1 nt-it])),3);
%     GG(:,:,it) = Numer(:,:,it)./Denom; 
% end

