% Note: savefast only works well with numerical arrays. see https://www.mathworks.com/matlabcentral/fileexchange/39721-save-mat-files-more-quickly


function saveRcvData_ULM(RcvData)
%     tic
    savepath = evalin('base', 'savepath');
    Trans = evalin('base', 'Trans');
    Psmall = evalin('base', 'Psmall');
    Psmall.bufferIndex = Psmall.bufferIndex + 1;
    assignin('base', 'Psmall', Psmall);
    paramNames = [num2str(round(Psmall.maxAngle)), '-', num2str(Psmall.na), '-', num2str(round(Psmall.frameRate)), '-', num2str(Psmall.numFramesPerBuffer), '-1-', num2str(Psmall.bufferIndex)];

%     filename = strcat(Trans.name, '-RcvData-', paramNames, '.mat');
    filename = strcat('RF-', paramNames, '.mat');

%     disp(strcat('Saving ', filename))

%     save([savepath, filename], 'RcvData', '-v7.3')
%     RcvData = RcvData(:, Trans.Connector, :);
%     whos RcvData
%     save([savepath, filename], 'RcvData', '-v7.3','-nocompression')
    savefast([savepath, filename], 'RcvData')


%     % new 1/15/25
%     fn = strcat("RcvData_buffer", num2str(Psmall.bufferIndex));
%     if Psmall.bufferIndex == 1
%         save([savepath, filename], fn, "-v6")
% 
%     else
%         save([savepath, filename], fn, "-v6", "-append")
%     end
%     toc
end