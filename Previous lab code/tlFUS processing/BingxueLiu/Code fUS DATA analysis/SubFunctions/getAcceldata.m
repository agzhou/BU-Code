function Acc = getAcceldata(Stim0)
% Stim0 is a vector of stimulation pattern;
% Stim0 = 0; resting state without stimulation;
% last modified by Bingxue Liu 03/04 2021;

        [accFileName,accFilePath]=uigetfile('H:\Accel data');
        accBas = readvars([accFilePath, accFileName], 'sheet',2);%'_unnamedTask<105>');
        if Stim0 == 0
            accStim = [];
        else
        accStim = readvars([accFilePath, accFileName], 'sheet',3);%'_unnamedTask<109>');
        end
        accRaw = [accBas; accStim];
        np = size(accRaw,1); % n points in all data
        nw  = 1000; % n points per window
        nbin = np/nw;
%         Stim0 = downsample(Stim0, size(Stim0,2)/nbin);
%        Stim0 = interp1(1:size(Stim0,2),Stim0,1:nbin);
        
       
        accRaw = sgolayfilt([0;diff(accRaw)],3,31); title('Raw Accel Data');
         figure; plot((1:np)/1000, accRaw);
        upthreshold = mean(accRaw)+10*std(accRaw);
        downthreshold = mean(accRaw)-10*std(accRaw);
        accRaw(find(accRaw>=upthreshold)) = upthreshold;
        accRaw(find(accRaw<=downthreshold)) = downthreshold;
     figure; plot((1:np)/1000, accRaw);
      hold on; plot((1:np)/1000,1.3*1e4*mean(diff(accRaw))*ones(np,1),'-.r'); title('Diff Accel Data');
      
       
        accBinVar = var(reshape(accRaw, [nw, nbin]),1);
        accSumVar = var(accRaw(:));
        ibin = find(accBinVar>accSumVar);
        figure; plot(1:nbin, accBinVar);
        hold on; plot(accSumVar*ones(nbin,1),'-.r');
        hold on; plot(ibin,accBinVar(ibin),'ro');
        hold on; area(1:nbin, Stim0*max(accBinVar).*ones(1,nbin),'FaceColor','r','EdgeColor','none'); alpha(0.1);
        xlabel('Time(Sec)');
        ylabel('Variance of Accel Data');
        Acc.accRaw = accRaw;
        Acc.accBinVar = accBinVar;
        Acc.accSumVar = accSumVar;
        Acc.ibin = ibin;
        Acc.nbin = nbin;
end
        




