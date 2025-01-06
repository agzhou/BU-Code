% Note: savefast only works well with numerical arrays. see https://www.mathworks.com/matlabcentral/fileexchange/39721-save-mat-files-more-quickly


% function saveRcvData(ReceiveData)
function saveRcvData_quick(RcvData)
%     tic
    savepath = evalin('base', 'savepath');
%     Trans = evalin('base', 'Trans');
%     Psmall = evalin('base', 'Psmall');
%     Psmall.supFrameIndex = Psmall.supFrameIndex + 1;
%     assignin('base', 'Psmall', Psmall);
%     paramNames = [num2str(Psmall.maxAngle), '-', num2str(Psmall.na), '-', num2str(Psmall.fps_target), '-', num2str(Psmall.numSubFrames), '-', num2str(Psmall.numSupFrames), '-', num2str(Psmall.supFrameIndex)];

    filename = strcat('RcvData-', datestr(now, 'mm-dd-yy-MM-SS'), '.mat');
%     disp(strcat('Saving ', filename))

%     save([savepath, filename], 'RcvData', '-v7.3')
%     RcvData = RcvData(:, Trans.Connector, :);
%     whos RcvData
    
    savefast([savepath, filename], 'RcvData')
%     toc
end