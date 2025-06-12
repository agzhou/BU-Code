clear all;
[FileName,FilePath]=uigetfile('F:\0320_2021_BL2\CP2\*.mat');
FileInfo = strsplit(FileName,'-');
newFilePath = 'E:\0423_BL3_dIQ\';
for i = 1:363
 thisFileName = strcat(strjoin(FileInfo(1:7),'-'),'-',num2str(i),'.mat');
 load([FilePath,thisFileName]);
% process data
 dIQ = diff(IQ,1,3);
%  newIQ = cat(3,dIQ(:,:,1),dIQ);
 IQ = dIQ;
 P.numCCframes = P.numCCframes-1;
 % save data;
 newthisFileName = strcat(strjoin(FileInfo(1:4),'-'),'-',num2str(P.numCCframes),'-',strjoin(FileInfo(6:7),'-'),'-',num2str(i),'.mat');
 save([newFilePath,newthisFileName],'IQ','P');
 disp(['save file ',newthisFileName])
 clear IQ P;
end

