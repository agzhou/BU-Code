% Note: savefast only works well with numerical arrays. see https://www.mathworks.com/matlabcentral/fileexchange/39721-save-mat-files-more-quickly


% function saveRcvData(ReceiveData)
function saveRcvData(RcvData)
    tic
    savepath = evalin('base', 'savepath');
%     Trans = evalin('base', 'Trans');
%     Psmall = evalin('base', 'Psmall');
    supFrameIndex = evalin('base', 'supFrameIndex');
%     Psmall.supFrameIndex = Psmall.supFrameIndex + 1;
    supFrameIndex = supFrameIndex + 1;
%     assignin('base', 'Psmall', Psmall);
    assignin('base', 'supFrameIndex', supFrameIndex);
%     paramNames = [num2str(Psmall.maxAngle), '-', num2str(Psmall.na), '-', num2str(Psmall.fps_target), '-', num2str(Psmall.numSubFrames), '-', num2str(Psmall.numSupFrames), '-', num2str(Psmall.supFrameIndex)];

%     filename = strcat(Trans.name, '-RcvData-', paramNames, '.mat');

    filename = strcat('RcvData-', num2str(supFrameIndex), '.mat');
%     disp(strcat('Saving ', filename))

%     save([savepath, filename], 'RcvData', '-v7.3')
    
    savefast([savepath, filename], 'RcvData')
    toc
end