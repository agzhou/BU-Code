% From previous lab members, modified a little

 function outputData = generateStimulus_variablerate(Absolute_Beginning_Stimulus_Delay_Seconds,Stimulus_Frequency_Hz...
    ,SingleStimulusPulseDurationSeconds,StimulateForHowManySeconds_InOneSequence, Total_One_Sequence_Duration_Seconds, daqrate)
    
    % Absolute_Beginning_Stimulus_Delay_Seconds = 1000;
    % Stimulus_Frequency_Hz = 3;
    % SingleStimulusPulseDurationSeconds = 1;
    % StimulateForHowManySeconds_InOneSequence = 5;
    % Total_One_Sequence_Duration_Seconds = 30;
    
    Absolute_Beginning_Stimulus_Delay_Seconds = Absolute_Beginning_Stimulus_Delay_Seconds./1000;
    SingleStimulusPulseDurationSeconds = SingleStimulusPulseDurationSeconds./1000;
 % Settings for CCD trigger
     Absolute_Beginning_Stimulus_Delay_Seconds_ch2  = 0;
     Stimulus_Frequency_Hz_ch2 = 15;
     SingleStimulusPulseDurationSeconds_ch2 = 0.05;
     StimulateForHowManySeconds_InOneSequence_ch2 = Total_One_Sequence_Duration_Seconds;
%      Total_One_Sequence_Duration_Seconds_ch2 = 30;
%      Number_Of_Sequence_Repetitions_ch2 = str2num(answer{6}) ;
%      Baseline_Delay_ch2 = 0;

    % set the output DAQ parameters...
    Output_DAQ_Frequency_Hz=daqrate;% sample rate is 10 kHz for output 
    Total_Number_Points = Total_One_Sequence_Duration_Seconds * Output_DAQ_Frequency_Hz;% * Number_Of_Sequence_Repetitions; % N

    % parameters for channel 2 - CCD trigger
    Number_Points_CCD_Trigger_Beginning_Delay = round(0*Output_DAQ_Frequency_Hz); % M - delay in triggering of the CCD - have a sense to be ZERO
     O = round(0.01*Output_DAQ_Frequency_Hz);% number of DAQ points in 10 ms duration - length of CCD trigger
    Number_Points_Absoulte_Beginning_ch2 = (Absolute_Beginning_Stimulus_Delay_Seconds_ch2*Output_DAQ_Frequency_Hz); % P
    Number_Points_For_Single_Stimulus_Period_ch2 = (Output_DAQ_Frequency_Hz/Stimulus_Frequency_Hz_ch2); % T
    Number_Points_For_Single_Stimulus_ON_ch2 = (SingleStimulusPulseDurationSeconds_ch2*Output_DAQ_Frequency_Hz); % R
    Number_Stimuluses_In_One_Sequence_ch2 = (StimulateForHowManySeconds_InOneSequence_ch2*Stimulus_Frequency_Hz_ch2); %k

    OutputChannel1 = zeros(Total_Number_Points,1);    
    OutputChannel2 = zeros(Total_Number_Points,1);


    % parameters for channel 1 - stimulation pattern
    Number_Points_Absolute_Beginning = (Absolute_Beginning_Stimulus_Delay_Seconds*Output_DAQ_Frequency_Hz); % P
    Number_Points_For_Single_Stimulus_Period = (Output_DAQ_Frequency_Hz/Stimulus_Frequency_Hz); % T
    Number_Points_For_Single_Stimulus_ON = (SingleStimulusPulseDurationSeconds*Output_DAQ_Frequency_Hz); % R
    Number_Stimuluses_In_One_Sequence = (StimulateForHowManySeconds_InOneSequence*Stimulus_Frequency_Hz); %k

    OutputChannel1 = zeros(Total_Number_Points,1);    
    OutputChannel2 = zeros(Total_Number_Points,1);
    
    % set CCD trigger channel values
%        OutputChannel1(Number_Points_CCD_Trigger_Beginning_Delay+1:Number_Points_CCD_Trigger_Beginning_Delay+O)= 5.0; % create 'O' duration square pulse, delayed for Number_Points_Stimulus_Beginning_Delay
%    OutputChannel1(Total_Number_Points) = 0.00001; % end the array with 0.00001 ???
%    pa1 = Number_Points_Absoulte_Beginning + 1; % skip absolute beginning delay        
%     for i = 1 : Number_Stimuluses_In_One_Sequence, % loop for each individual stimulus
%         OutputChannel1(pa1:pa1+Number_Points_For_Single_Stimulus_ON-1) = 5.0; % square pulse for each individual stimulus
%         pa1 = pa1 + Number_Points_For_Single_Stimulus_Period; % skip to the next individual stimulus
%     end;

    
    %Set one CCD sequence channel values
   OutputChannel1(Number_Points_CCD_Trigger_Beginning_Delay+1:Number_Points_CCD_Trigger_Beginning_Delay+O)= 5.0; % create 'O' duration square pulse, delayed for Number_Points_Stimulus_Beginning_Delay
   OutputChannel1(Total_Number_Points) = 0.00001; % end the array with 0.00001 ???
   pa1 = Number_Points_Absoulte_Beginning_ch2 + 1; % skip absolute beginning delay        
    for i = 1 : Number_Stimuluses_In_One_Sequence_ch2 % loop for each individual stimulus
        OutputChannel1(pa1:pa1+Number_Points_For_Single_Stimulus_ON_ch2-1) = 5.0; % square pulse for each individual stimulus
        pa1 = pa1 + Number_Points_For_Single_Stimulus_Period_ch2; % skip to the next individual stimulus
    end
    % set one stimulation sequence channel values   
    pa1 = Number_Points_Absolute_Beginning + 1; % skip absolute beginning delay        
    for i = 1 : Number_Stimuluses_In_One_Sequence % loop for each individual stimulus
        OutputChannel2(pa1:pa1+Number_Points_For_Single_Stimulus_ON-1) = 5.0; % square pulse for each individual stimulus
        pa1 = pa1 + Number_Points_For_Single_Stimulus_Period; % skip to the next individual stimulus
    end
    OutputChannel1(Total_Number_Points)=0.00001;
    OutputChannel2(Total_Number_Points)=0.00001; % end the array again with 0.00001 ???
    outputData(:,1) = OutputChannel2;
    outputData(:,2) = OutputChannel1;
    end