% Note: savefast only works well with numerical arrays. see https://www.mathworks.com/matlabcentral/fileexchange/39721-save-mat-files-more-quickly


function saveRcvData_ULM(RcvData)
%     tic
    savepath = evalin('base', 'savepath');
    Trans = evalin('base', 'Trans');
    Psmall = evalin('base', 'Psmall');
    Psmall.bufferIndex = Psmall.bufferIndex + 1;
    assignin('base', 'Psmall', Psmall);
    paramNames = [num2str(Psmall.maxAngle), '-', num2str(Psmall.na), '-', num2str(Psmall.PRF), '-', num2str(Psmall.frameRate), '-', num2str(Psmall.numFramesPerBuffer), '-', num2str(Psmall.bufferIndex)];

    filename = strcat(Trans.name, '-RcvData-', paramNames, '.mat');
%     disp(strcat('Saving ', filename))

%     save([savepath, filename], 'RcvData', '-v7.3')
%     RcvData = RcvData(:, Trans.Connector, :);
%     whos RcvData
    
    savefast([savepath, filename], 'RcvData')
%     toc
end