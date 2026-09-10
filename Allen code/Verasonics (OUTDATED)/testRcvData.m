function testRcvData(ReceiveData)

%     savepath = evalin('base', 'savepath');
%     Trans = evalin('base', 'Trans');
%     RcvData = evalin('base', 'RcvData');
    RcvData = ReceiveData;

    disp(all(all(RcvData == 0)))

%     date_time = char(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss.SSS'));
% %     date_time(date_time == '-') = '_';
% %     date_time(date_time == ':') = '_';
% %     date_time(date_time == ' ') = '_';
%     date_time(date_time == '-' | date_time == ':' | date_time == ' ' | date_time == '.') = '_';
% %     filename = strcat('RC15gV_RcvData_', date_time, '.mat');
%     filename = strcat(Trans.name, '_RcvData_', date_time, '.mat');
% %     disp(strcat('Saving ', filename))
% 
%     save([savepath, filename], 'RcvData', '-v7.3')

end