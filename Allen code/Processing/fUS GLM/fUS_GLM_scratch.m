

% trange - defines the range for the block average [tPre tPost dt]. If dt
%           defined, time series are interpolated prior to 

function [] = fUS_GLM_scratch(data, t, stim, Aaux, tIncAuto, trange, glmSolveMethod, idxBasis, paramsBasis, driftOrder)

    % I don't want to use the DataClass and snirf stuff

    %%%% Get characteristics of the timing
    dt = t(2) - t(1); % Time step %%%%%%%%%%%%%
    fq = 1/dt;        % Sampling/measurement frequency
    nPre = round(trange(1)/dt); % Which index for the block average Pre point
    nPost = round(trange(2)/dt);% Which index for the block average Post point
    % nTpts = size(y,1); % # of time points
    nTpts = length(t); % # of time points
    tHRF = (1*nPre*dt:dt:nPost*dt)';
    ntHRF = length(tHRF); % # of points in the HRF
    nT = length(t); % Also # of time points?

    %%%% Update stimulus vectors
    % Note: it's assumed that the stimulus vector has the same length and
    % timesteps of the data vector.
    % stim represents the onsets of the active stimulus periods.
    stimStates = stim;
    stimStates(stimStates > 0) = stimStates(stimStates > 0) ./ stimStates(stimStates > 0);
    stimStates = logical(stimStates); % stimStates: a logical version of the stim vector. True = stim on, false = stim off
    stimAmps = stim; % stimAmps: a quantitative version of the stim vector. The amplitude of the stimulus is preserved.
    
    
    
    % for 


    %%%%%%%%%%%%%%%%
    % Prune good stim, generate onset matrix
    %%%%%%%%%%%%%%%%
    % Get only indices of conditions with any stimStates that are 1
    % lstCond = find(sum(stimStates == 1, 1) > 0);
    lstCond = find(stimStates == true); % Indices at the active stim onsets
    % nCond = length(lstCond); % # of stim onsets
    % nTrials = zeros(nCond, 1); % # of trials??......
    nCond = 1; % Only one condition
    nTrials = length(lstCond);

    onset = zeros(nT, nCond);
    avg_pulses = {};
    for iCond = 1:nCond
        % lstT = find(stimStates(:, lstCond(iCond)) == 1);  % Indices of stims enabled (== 1)
        % lstp = find((lstT+nPre) >= 1 & (lstT+nPost) <= nTpts);  % Indices of stims not clipped by signal
        lstT = find(stimStates == true);  % Indices of stims enabled (== 1)
        lstp = find((lstT+nPre) >= 1 & (lstT+nPost) <= nTpts);  % Indices of stims not clipped by signal
        lst = lstT(lstp); % Final list of stim onset indices to use
        % nTrials(iCond) = length(lst);

        % Generate basis boxcars of stim amplitude and duration
        starts = lst+nPre; % Get the start indices of each trial's "stim" --> include the Pre offset
        % if ~isempty(stim(lstCond(iCond)))
            durations = stim(lstCond(iCond)).data(:, 2);
            amplitudes = stim(lstCond(iCond)).data(:, 3);
            avg_pulses{iCond} = ones(round(mean(durations) / dt), 1); %#ok<AGROW>
            for i = 1:length(starts)
                if idxBasis == 1  % Gaussian has no duration T (yet)
                   pulse_duration = 1; 
                else
                   pulse_duration = round(durations(i) / dt); 
                end
                pulse = (amplitudes(i) / pulse_duration) * ones(pulse_duration, 1);
                onset(starts(i):starts(i) + pulse_duration - 1, iCond) = onset(starts(i):starts(i) + pulse_duration - 1, iCond) + pulse;
            end
        % end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct the basis functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch idxBasis
        case 1
            % Gaussians
            gms = paramsBasis(1);
            gstd = paramsBasis(2);
            
            nB = floor((trange(2)-trange(1)) / gms) - 1;
            tbasis = zeros(ntHRF,nB);
            for b=1:nB
                tbasis(:,b) = exp(-(tHRF-(trange(1)+b*gms)).^2/(2*gstd.^2));
                tbasis(:,b) = tbasis(:,b)./max(tbasis(:,b));
            end
        
        case 2
            % Modified Gamma
            if length(paramsBasis)==2
                nConc = 1;
            else
                nConc = 2;
            end
            
            nB = 1;
            tbasis = zeros(ntHRF,nB,nConc);
            for iConc = 1:nConc
                tau = paramsBasis((iConc-1)*2+1);
                sigma = paramsBasis((iConc-1)*2+2);
                
                tbasis(:,1,iConc) = (exp(1)*(tHRF-tau).^2/sigma^2) .* exp( -(tHRF-tau).^2/sigma^2 );
                lstNeg = find(tHRF<0);
                tbasis(lstNeg,1,iConc) = 0;
                
                if tHRF(1)<tau
                    tbasis(1:round((tau-tHRF(1))/dt),1,iConc) = 0;
                end
                
            end
        
        case 3
            % Modified Gamma and Derivative
            if length(paramsBasis)==2
                nConc = 1;
            else
                nConc = 2;
            end
            
            nB = 2;
            tbasis=zeros(ntHRF,nB,nConc);
            for iConc = 1:nConc
                tau = paramsBasis((iConc-1)*2+1);
                sigma = paramsBasis((iConc-1)*2+2);
                
                tbasis(:,1,iConc) = (exp(1)*(tHRF-tau).^2/sigma^2) .* exp( -(tHRF-tau).^2/sigma^2 );
                tbasis(:,2,iConc) = 2*exp(1)*( (tHRF-tau)/sigma^2 - (tHRF-tau).^3/sigma^4 ) .* exp( -(tHRF-tau).^2/sigma^2 );
                
                if tHRF(1)<tau
                    tbasis(1:round((tau-tHRF(1))/dt),1:2,iConc) = 0;
                end
                
            end
        
        case 4
            % AFNI Gamma function
            if length(paramsBasis)==2
                nConc = 1;
            else
                nConc = 2;
            end
            
            nB=1;
            tbasis=zeros(ntHRF,nB,nConc);
            for iConc = 1:nConc
                
                p = paramsBasis((iConc-1)*2+1);
                q = paramsBasis((iConc-1)*2+2);
                
                tbasis(:,1,iConc) = (tHRF/(p*q)).^p.* exp(p-tHRF/q);
                
            end
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct design matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dA=zeros(nT,nB*nCond,2);
    for iConc = 1:2
        iC = 0;
        for iCond=1:nCond
            for b=1:nB
                iC = iC + 1;
                if size(tbasis,3)==1
                    clmn = conv(onset(:,iCond),tbasis(:,b));
                else
                    clmn = conv(onset(:,iCond),tbasis(:,b,iConc));
                end
                clmn = clmn(1:nT);
                dA(:,iC,iConc) = clmn;
                beta_label{b + (iCond-1)*nB} = ['Cond' num2str(iCond)];
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Expand design matrix nth order polynomial for drift correction
    % rescale polynomial to avoid bad conditionning
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xDrift = ones(nT,driftOrder);
    for ii=2:(driftOrder+1)
        xDrift(:,ii) = ([1:nT]').^(ii-1);
        xDrift(:,ii) = xDrift(:,ii) / xDrift(end,ii);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Expand design matrix with Aaux
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nAux = size(Aaux,2);
    if flagNuisanceRMethod == 3
        nAux = 0;
    end

end