 clear all; 

[FileName,FilePath]=uigetfile('G:\0622_BL3_sedated\');
fileInfo=strsplit(FileName(1:end-4),'-');
myFile=matfile([FilePath,FileName]);
P=myFile.P;
prompt={'Start Repeat', 'Number of Repeats','Time Interval [s]'};
name='File info'; 
defaultvalue={'31','32','1',num2str(1/P.PDIFR-P.numCCframes/P.CCFR)};
numinput=inputdlg(prompt,name, 1, defaultvalue);
startRpt=str2num(numinput{1});
nRpt=str2num(numinput{2});          % number of repeat for each coronal plane
tIntV=str2num(numinput{3});        % time interval [s] between two IQ datasets

indSkipped=1; k = 1; 
mIQ = single(zeros(length(P.zCoor),length(P.xCoor),nRpt*P.numCCframes));
for iRpt=startRpt:startRpt+nRpt-1
 iFileInfo=fileInfo;
 iFileInfo{8}=num2str(iRpt);
 iFileName=[strjoin(iFileInfo,'-'),'.mat'];
            %     load([FilePath,iFileName]);
 if exist([FilePath,iFileName],'file')
  myFile=matfile([FilePath,iFileName]);
  IQ = myFile.IQ;
  disp([iFileName,' was loaded!'])
  else
   disp([iFileName, ' skipped!'])
  SkipFile(indSkipped)=iRpt;
   indSkipped=indSkipped+1;  
 end
 IQ = single(IQ);
  mIQ(:,:,((k-1)*P.numCCframes+1):k*P.numCCframes) = IQ;
  k = k+1;
end
% P.numCCframes = nRpt*P.numCCframes;
% P.tIntV = tIntV;
% iFileInfo{5} = num2str(P.numCCframes);
% save([FilePath,strjoin(iFileInfo,'-'),'.mat'],'-v7.3','mIQ','P');

sumIQ = squeeze(sum(sum(mIQ,1),2));
tcoor = linspace(1,150,150*200);
figure;subplot(411); plot(tcoor,abs(sumIQ)); xlabel('s');title('abs(sum(IQ))');
subplot(412); plot(tcoor,real(sumIQ)); xlabel('s');title('real(sum(IQ))');
subplot(413); plot(tcoor,imag(sumIQ)); xlabel('s');title('imag(sum(IQ))');
subplot(414); plot(tcoor,angle(sumIQ));xlabel('s');title('phase(sum(IQ))');

        
        