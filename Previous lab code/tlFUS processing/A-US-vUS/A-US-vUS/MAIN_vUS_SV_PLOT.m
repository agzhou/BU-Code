%% load and plot Power Doppler results
%% for multiple files PDI process
clear all;
defaultpath='Z:\US-DATA\PROJ-D-vUS\Phantom-Validation-20181016-FlowRBC\DATA-FlowRBC-20181016-AngledX\';
[FileName,FilePath]=uigetfile(defaultpath);  % read data of a small part of the brain cortex (IQR matrix)  
fileInfo=strsplit(FileName(1:end-4),'-');
startCP0=regexp(fileInfo{7},'\d*','Match');
cpName=regexp(fileInfo{7},'\D*','Match');
if isempty(startCP0{1})
    startCP1='0';
else
    startCP1=startCP0{1};
end
% data processing parameters
prompt={'Start file (CP)','Number of files (CPs)','Start Repeat', 'Number of Repeats','nRfn','Angle(0),Trans(1)'};
name='File info';
defaultvalue={num2str(startCP1), '11',fileInfo{8},'1','1','0'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
startCP=str2num(numinput{1});
nCP=str2num(numinput{2});          % number of coronal planes
startRpt=str2num(numinput{3});
nRpt=str2num(numinput{4});          % number of repeat for each coronal plane
nRfn=str2num(numinput{5});          % image refine scale
dataType=str2num(numinput{6}); 
[VzCmap,VzCmapDn, VzCmapUp, pdiCmapUp, PhtmCmap]=Colormaps_fUS;
%% load file and plot
vAct=[1:2:15 20 25 30];
for iCP=startCP:startCP+nCP-1
    for iRpt=startRpt:startRpt+nRpt-1
        iFileInfo=fileInfo;
        if nCP>1
            iFileInfo{7}=[cpName{1},num2str(iCP)];
        end              
        iFileInfo{8}=num2str(iRpt);
        iFileName=[strjoin(iFileInfo,'-'),'.mat'];
%         load([FilePath,iFileName]);
        myFile=matfile([FilePath,iFileName]);
        Vx=myFile.Vx;
        Vz=myFile.Vz;
        Vcz=myFile.Vcz;
        R=myFile.R;
        Mf=myFile.Mf;
        P=myFile.P;
        PDI=(myFile.PDI)/2e15;
        [nz,nx,nPDI]=size(Vz);
        %% Vx weighted with R
        Rmsk=0.01*ones(size(R));
        Rmsk(R>0.2)=R(R>0.2);
        Rmsk(Rmsk>0.6)=1;
%         VxR=Vx(:,:,1:2).*Rmsk(:,:,1:2);
        
        MfR=Mf(:,:,2);
        MfRmsk=0.01*ones(size(MfR));
        MfRmsk(MfR>0.1)=MfR(MfR>0.1);
        MfRmsk(MfRmsk>0.4)=1;
        VxR=Vx.*Rmsk;
        
        Vx=Vx.*Rmsk; Vz=Vz.*Rmsk;
        %% V total
        V=sign(Vz).*sqrt(Vx.^2+Vz.^2);
        %% Refine image and plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [nzPDI,nxPDI,~]=size(PDI);
        PDIn=zeros(size(PDI));
        PDIrfn=zeros(nz*nRfn,nx*nRfn,2);
        Vrfn=zeros(nz*nRfn,nx*nRfn,2);
        Vzrfn=zeros(nz*nRfn,nx*nRfn,2);
        
        if dataType==0
            %% angled flow phantom velocity in the tube
            [CX CY C]=improfile(V,[8 247], [13 115]);
            for iP=1:240
                Vroi(:,:,iP)=abs(V(floor(CY(iP))+[-3:3],floor(CX(iP))));
                Vzroi(:,:,iP)=abs(Vz(floor(CY(iP))+[-3:3],floor(CX(iP))));
                Vxroi(:,:,iP)=abs(Vx(floor(CY(iP))+[-3:3],floor(CX(iP))));
                Vczroi(:,:,iP)=abs(Vcz(floor(CY(iP))+[-3:3],floor(CX(iP))));
                PDIroi(:,:,iP)=abs(PDI(floor(CY(iP))+[-3:3],floor(CX(iP)),3));
            end
            Vmean=mean(Vroi(:));
            Vstd=std(Vroi(:));
            Vzmean=mean(Vzroi(:));
            Vzstd=std(Vzroi(:));
            Vxmean=mean(Vxroi(:));
            Vxstd=std(Vxroi(:));
            Vczmean=mean(Vczroi(:));
            Vczstd=std(Vczroi(:));
            PDImean=mean(PDIroi(:));
            PDIstd=std(PDIroi(:));
        else
            %% Transverse phantom velocity in the tube
            Vroi=abs(V(54:57,:));
            Vmean=mean(Vroi(:));
            Vstd=std(Vroi(:));
            Vzroi=abs(Vz(54:57,:));
            Vzmean=mean(Vzroi(:));
            Vzstd=std(Vzroi(:));
            Vxroi=abs(Vx(54:57,:));
            Vxmean=mean(Vxroi(:));
            Vxstd=std(Vxroi(:));
            Vczroi=abs(Vcz(54:57,:));
            Vczmean=mean(Vczroi(:));
            Vczstd=std(Vczroi(:));
            PDIroi=abs(PDI(54:57,:,3));
            PDImean=mean(PDIroi(:));
            PDIstd=std(PDIroi(:));
        end
        Vall(:,iCP,iRpt,1)=[Vmean,Vstd]; % save all mean V and std V
        Vall(:,iCP,iRpt,2)=[Vzmean,Vzstd]; % save all mean Vz and std Vz
        Vall(:,iCP,iRpt,3)=[Vxmean,Vxstd]; % save all mean Vx and std Vx
        Vall(:,iCP,iRpt,4)=[Vczmean,Vczstd]; % save all mean Vcz and std Vcz
        Vall(:,iCP,iRpt,5)=[PDImean,PDIstd]; % save all mean PDI and std PDI
        %% figure plot
        Fig=figure;
        set(Fig,'Position',[400 400 900 400])
        subplot(2,2,1);imagesc(P.xCoor, P.zCoor,abs(V));
        colormap(PhtmCmap);
        caxis([0 30])
        title(['Vset=',num2str(vAct(iCP)),' mm/s, Vrec=',num2str(Vmean),'+/-', num2str(Vstd),' mm/s'])
        xlabel('X [mm]')
        ylabel('Y [mm]')
        colorbar
        axis equal tight
        
        
        hAxes1=subplot(2,2,2);
        imagesc(P.xCoor, P.zCoor,abs(PDI(:,:,3)));
        colormap(hAxes1,hot);
        caxis([0 1])
        title(['Vset=',num2str(vAct(iCP)),' mm/s, PDI=',num2str(PDImean),'+/-', num2str(PDIstd)])
        xlabel('X [mm]')
        ylabel('Y [mm]')
        colorbar
        axis equal tight
        
        subplot(2,2,3);imagesc(P.xCoor, P.zCoor,abs(Vz));
        colormap(PhtmCmap);
        caxis([0 30])
        title(['Vset=',num2str(vAct(iCP)),' mm/s, Vzrec=',num2str(Vzmean),'+/-', num2str(Vzstd),' mm/s'])
        xlabel('X [mm]')
        ylabel('Y [mm]')
        colorbar
        axis equal tight
        
        subplot(2,2,4);imagesc(P.xCoor, P.zCoor,abs(Vcz));
        colormap(PhtmCmap);
        caxis([0 30])
        title(['Vset=',num2str(vAct(iCP)),' mm/s, Vcz=',num2str(Vczmean),'+/-', num2str(Vczstd),' mm/s'])
        xlabel('X [mm]')
        ylabel('Y [mm]')
        colorbar
        axis equal tight
        
        saveas(gcf,[FilePath,'V',iFileName(4:end-4),'.png'],'png');
        saveas(gcf,[FilePath,'V',iFileName(4:end-4),'.fig'],'fig');
        
    end
end
save([FilePath,'statisticsVnPDI.mat'],'Vall');
%% statistic result
Vset=[1:2:15 20 25 30];
P=polyfit(Vset,squeeze(mean(Vall(1,:,iRpt,1),3)),1);
Y=polyval(P,Vset);
figure;
yyaxis left;
errorbar(Vset,squeeze(mean(Vall(1,:,iRpt,1),3)),squeeze(mean(Vall(2,:,iRpt,1),3)),'k.');hold on, plot(Vset,Y,'k');
text(15,13,['v=',num2str(P(1),'%2.2f'),'v_s_e_t',num2str(P(2),'%2.2f')],'FontSize',10)
axis equal
ylim([0 35])
xlim([0 35])

P=polyfit(Vset,squeeze(mean(Vall(1,:,iRpt,3),3)),1);
Y=polyval(P,Vset);
hold on, errorbar(Vset,squeeze(mean(Vall(1,:,iRpt,3),3)),squeeze(mean(Vall(2,:,iRpt,3),3)),'g.');hold on, plot(Vset,Y,'g');
text(15,13,['v_x=',num2str(P(1),'%2.2f'),'v_s_e_t',num2str(P(2),'%2.2f')],'FontSize',10)
axis equal
ylim([0 35])
xlim([0 35])

P=polyfit(Vset,squeeze(mean(Vall(1,:,iRpt,2),3)),1);
Y=polyval(P,Vset);
hold on, errorbar(Vset,squeeze(mean(Vall(1,:,iRpt,2),3)),squeeze(mean(Vall(2,:,iRpt,2),3)),'b.');hold on, plot(Vset,Y,'b');
text(15,13,['v_z=',num2str(P(1),'%2.2f'),'v_s_e_t',num2str(P(2),'%2.2f')],'FontSize',10)
axis equal
ylim([0 35])
xlim([0 35])

P=polyfit(Vset,squeeze(mean(Vall(1,:,iRpt,4),3)),1);
Y=polyval(P,Vset);
hold on, errorbar(Vset,squeeze(mean(Vall(1,:,iRpt,4),3)),squeeze(mean(Vall(2,:,iRpt,4),3)),'.');hold on, plot(Vset,Y);
text(15,13,['v_z_c_D=',num2str(P(1),'%2.2f'),'v_s_e_t',num2str(P(2),'%2.2f')],'FontSize',10)
axis equal
ylim([0 35])
xlim([0 35])

hold on, plot(Vset,squeeze(mean(Vall(1,:,iRpt,1),3)),'k.');
hold on, plot(Vset,squeeze(mean(Vall(1,:,iRpt,2),3)),'b.');
hold on, plot(Vset,squeeze(mean(Vall(1,:,iRpt,3),3)),'g.');
hold on, plot(Vset,squeeze(mean(Vall(1,:,iRpt,4),3)),'b.');
yyaxis right;
errorbar(Vset,squeeze(mean(Vall(1,:,iRpt,5),3)),squeeze(mean(Vall(2,:,iRpt,5),3)),'r.');
plot(Vset,squeeze(mean(Vall(1,:,iRpt,5),3)),'r.-');
% figure;
% errorbar(Vset,squeeze(mean(Vall(1,:,5,1),3)),squeeze(mean(Vall(2,:,5,1),3)));
% hold on, errorbar(Vset,squeeze(mean(Vall(1,:,5,2),3)),squeeze(mean(Vall(2,:,5,2),3)));
% hold on, errorbar(Vset,squeeze(mean(Vall(1,:,5,3),3)),squeeze(mean(Vall(2,:,5,3),3)));
% hold on, errorbar(Vset,squeeze(mean(Vall(1,:,5,4),3)),squeeze(mean(Vall(2,:,5,4),3)));
% hold on, errorbar(Vset,squeeze(mean(Vall(1,:,5,5),3)),squeeze(mean(Vall(2,:,5,5),3)));
% axis equal
% xlim([0 35])
% ylim([0 35])
% grid on
