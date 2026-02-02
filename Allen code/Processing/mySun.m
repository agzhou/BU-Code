
% cd("C:\Users\BOAS-US\Documents\Allen\GitHub\BU-Code\Allen code\Processing")

IQData = squeeze(complex(IData,QData));
chk2 = permute(squeeze(max(IQData,[],1)),[2,1,3]);
%%

p.szAcq = P.Receive(1).endSample;
p.ConnMap = P.Trans.Connector;
p.na = P.na*2;
p.nFrames = P.numFramesPerBuffer;
p.numEl = P.numElements*2;

RFData = zeros(p.szAcq,p.numEl,p.na,p.nFrames,'int16');

for i = 1:p.nFrames
    for j = 1:p.na
        RFData(:,:,j,i) = RcvData(P.Receive(j).startSample:P.Receive(j).endSample,p.ConnMap,i);
    end
end




%% Plot RF XCorrs

cd("M:\Ultrasound data from 01-15-2026 to\01-15-2026 AZ01 FC RC15gV continuous run 2");
RFName = "RF-5-11-400-400-1-";
xCorrSignal = zeros(p.nFrames,10);
for i = 1:10
    name = RFName + num2str(i) + ".mat";

    RFData2 = load(name,'RcvData').('RcvData');
    RFData2 = RFData2(:,p.ConnMap,:);
    xCorrSignal(:,i) = calcRFXC(RFData2);
end

%%
figure
plot(xCorrSignal)



%% test

timeTags = readTimeTags(RcvData);

timeTags = timeTags - timeTags(1);



%% Helper functions
function [timeTags] = readTimeTags(RFData)
    
    timeTags = zeros(size(RFData,3),1);
    for frmCount = 1:size(RFData,3)
    
        tmp2 = [0;0];
        if (any(RFData(:,:,frmCount),'all'))
            tmp = RFData(1:2,:,frmCount);
            tmp2(1,1) = tmp(1,find(tmp(1,:),1));
            tmp2(2,1) = tmp(2,find(tmp(2,:),1));
            timeStamp = getTimeStamp(double(tmp2));
        end
        

        


%         timeStamp = getTimeStamp(double(RFData(1:2,1,frmCount)));
        % the 32 bit time tag counter increments every 25 usec, so we have to scale
        % by 25 * 1e-6 to convert to a value in seconds

        timeTags(frmCount) = timeStamp/4e4;
    end
end

function [tStamp] = getTimeStamp(W)

    % get time tag from first two samples
    % time tag is 32 bit unsigned interger value, with 16 LS bits in sample 1
    % and 16 MS bits in sample 2.  Note RDatain is in signed INT16 format so must
    % convert to double in unsigned format before scaling and adding
    for i=1:2
        if W(i) < 0
            % translate 2's complement negative values to their unsigned integer
            % equivalents
            W(i) = W(i) + 65536;
        end
    end
    tStamp = W(1) + 65536 * W(2);

end