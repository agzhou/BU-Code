% calculate 3D g1, return nz_nx_ny_ntau
% nt: number of time points for ACF calculation, usuall, nt=nxRpt
function [GG, Numer] = IQ2g1(IQ, Start_T, nt, ntau)
%% constant total number of samples, nt
[nz, nx, nxRpt] = size(IQ) ;
if nt>nxRpt-ntau
    nt=nxRpt-ntau-Start_T+1;
%     disp(['Warning: nt is larger than nxRpt-ntau, and is modified to be nxRpt-ntau=',num2str(nt),'!']);
end
%% calculate g1
% temp_deno=zeros(nz,nx,nt);
% temp_numer=zeros(nz,nx,ntau,nt);
% temp_deno(:,:,:)=(conj(IQ(:,:,Start_T:Start_T-1+nt))).*(IQ(:,:,Start_T:Start_T-1+nt));
% for itau = 1:ntau
%     temp_numer(:,:,itau,:)=(conj(IQ(:,:,Start_T:Start_T-1+nt))).*(IQ(:,:,itau+Start_T:itau+Start_T-1+nt));
% end
% Denom=repmat(mean(temp_deno,3),[1,1,ntau]); % calculate the denominator
% Numer=(mean(temp_numer,4));
% GG = Numer./Denom;
%% calculate g1 pixel by pixel
GG=zeros(nz,nx,ntau);
for iz=1:nz
    temp_deno=zeros(1,nx,nt);
    temp_numer=zeros(1,nx,ntau,nt);
    temp_deno=(conj(IQ(iz,:,Start_T:Start_T-1+nt))).*(IQ(iz,:,Start_T:Start_T-1+nt));
    for itau = 1:ntau
        temp_numer(1,:,itau,:)=(conj(IQ(iz,:,Start_T:Start_T-1+nt)).*(IQ(iz,:,itau+Start_T:itau+Start_T-1+nt)));
    end
    Denom=repmat(mean(temp_deno,3),[1,1,ntau]); % calculate the denominator
    Numer=(mean(temp_numer,4));
    GG(iz,:,:) = Numer./Denom;
%     figure;plot(abs(squeeze(Numer)));
end
% %% varying total number of samples, nt-itau
% ny=1;
% [nz, nx, nxRpt] = size(IQ) ;
% %%%% Jonghwan ACF mean((yi*)*(y(i+tau))/mean(yi**yi);
% temp=zeros(nz,nx,nt);
% Numer=zeros(nz,nx,ntau);
% for itt=Start_T:Start_T-1+nt
%     temp(:,:,itt-Start_T+1)=squeeze(IQ(:,:,itt))'.'.*squeeze(IQ(:,:,itt));
% end
% Denom=mean(temp,3); % calculate the denomenator
% for it = 1:ntau
%     %% calcualte numerator 
%     temp=zeros(nz,nx,nt);
%         for itt=Start_T:Start_T-1+nt-it
%             
%             temp(:,:,itt-Start_T+1)=squeeze(IQ(:,:,itt)'.').*squeeze(IQ(:,:,it+itt));
%         end
%     Numer(:,:,it)=mean(temp,3);
%     GG(:,:,it) = Numer(:,:,it)./Denom; 
% end

