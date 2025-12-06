function HRF = calROIHRF(cData,trial, BW)
nRest = trial.nRest;
nBase = trial.nBase;
nStim = trial.nStim;
nTrial = trial.n;
nRecover = trial.nRecover;
nlength = nStim+nRecover;
startRpt = 1;
nRpt = length(trial.stim);

if BW == 0
h=msgbox(['Plese select the responding ROI...']);
figure(2)
pause(0.5)
delete(h)
[xROI, zROI]=ginput(6);
xROI=floor(xROI);
zROI=floor(zROI);
BW=[];
BW(:,:,1)=roipoly(cData(:,:,1),xROI,zROI);
end
    %% Value in the selected ROI
    % baseline pixel index for the ROI t course calculation
    DataBase=BW(:,:).*mean(cData(:,:,1:trial.nBase),3);
    Locs=abs(DataBase)~=0;
    for iRpt=startRpt:startRpt+nRpt-1
        Data=BW(:,:).*cData(:,:,iRpt);
        roi=Data(Locs>0);
        mData(iRpt)=mean(roi(:));
        stdData(iRpt)=std(roi(:));
        semData(iRpt)=std(roi(:))/sqrt(numel(roi));
    end    

%     %% time course plot for each ROI
%     xCoor=[startRpt:startRpt+nRpt-1];%*handles.tIntV;
%     
    %% trial averaged time course
    xCoor=[1:nlength];%*handles.tIntV;
    stimTrial=trial.stim(nBase-nRest+1:nBase-nRest+nlength);
    
    Trial=reshape(mData(nBase-nRest+1:nRpt-nRest),[nlength,nTrial]);
    mTrial(1,:)=median(Trial,2)'; % mean
    stdTrial(1,:)=std(Trial,1,2)';
    semTrial(1,:)=(std(Trial,1,2)./sqrt(nTrial))';
    RatioTrial=(Trial-repmat(mean(Trial(1:nRest,:),1),[size(Trial,1),1]))./abs(mean(Trial(1:nRest,:),1));% mean);% 5
    mRatioTrial(1,:)=median(RatioTrial,2)'; % mean
    stdRatioTrial(1,:)=std(RatioTrial,1,2)';
    semRatioTrial(1,:)=(std(RatioTrial,1,2)./sqrt(nTrial))';

%% obtain response function
mResponse=mRatioTrial*100;
HRFmov=movmean(mResponse,5);% averaged ratio trial after movemean 
lb = [1     0.1     0.1     xCoor(5)]; % [a, b, c, t0]
ub = [200   10      2     xCoor(8)]; 
fun = fittype('a*(t-t0).^2.*exp(-b*abs(t-t0).^c).*normcdf(t,t0,0.006)','indep','t');
Fopts = fitoptions(fun);
Fopts.Lower = lb;
Fopts.Upper = ub;
Fopts.Display='off';
Fopts.TolFun=1e-6;
Fopts.Robust='LAR';
Fopts.StartPoint=[150, 2, 1, xCoor(5)];
est = fit(xCoor',mResponse',fun,Fopts);
HRF=fun(est.a,est.b,est.c,est.t0,xCoor);
R_HRF=(1-sum(abs(mResponse-HRF).^2)./sum(abs((mResponse)-mean(mResponse)).^2));
R_HRFmov=(1-sum(abs(mResponse-HRFmov).^2)./sum(abs((mResponse)-mean(mResponse)).^2));
COLOR=[1 0 0;
    0.08 0.17 0.55
    0.31 0.31 0.31];
tCoor=xCoor-xCoor(5);
Fig=figure;
set(Fig,'Position',[800 800 300 200])
    Color.Shade=COLOR(1,:);
    Color.ShadeAlpha=0.3;
    Color.Line=COLOR(1,:);
    %% avearged relative change
    hold on;     ShadedErrorbar(mRatioTrial(1,:)*100,semRatioTrial(1,:)*100,tCoor,Color); 
    title('Trial averaged response')
    ylabel('Relative Change [%]')
    xlabel('t [s]')
% hold on, plot(tCoor, stimTrial*(max(abs(mRatioTrial(:)))*100-100)+100,'b')
hold on, plot(tCoor, stimTrial*(max(mRatioTrial(:)))*100,'b');
set(gca, 'YGrid', 'on', 'XGrid', 'off')
hold on, plot(tCoor,HRF,'g','LineWidth',2)
title(['R=',num2str(R_HRF)])
