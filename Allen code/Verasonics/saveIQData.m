% Note: savefast only works well with numerical arrays. see https://www.mathworks.com/matlabcentral/fileexchange/39721-save-mat-files-more-quickly

%%%% NOTE: THIS FILE DOES NOT RUN BECAUSE runAcq WILL CALL FROM THE
%%%% VERASONICS FOLDER %%%%

% function saveRcvData(ReceiveData)
function saveIQData(IQData)
%     tic
    disp('test')
%     disp(size(IQData))

%     savepath = evalin('base', 'savepath');
%     Trans = evalin('base', 'Trans');
%     Psmall = evalin('base', 'Psmall');
%     Psmall.supFrameIndex = Psmall.supFrameIndex + 1;
%     assignin('base', 'Psmall', Psmall);
%     paramNames = [num2str(Psmall.maxAngle), '-', num2str(Psmall.na), '-', num2str(Psmall.fps_target), '-', num2str(Psmall.numSubFrames), '-', num2str(Psmall.numSupFrames), '-', num2str(Psmall.supFrameIndex)];
% 
%     filename = strcat(Trans.name, '-RcvData-', paramNames, '.mat');
% %     disp(strcat('Saving ', filename))
% 
% %     save([savepath, filename], 'RcvData', '-v7.3')
% %     RcvData = RcvData(:, Trans.Connector, :);
% %     whos RcvData
%     
%     savefast([savepath, filename], 'RcvData')

%     toc
end