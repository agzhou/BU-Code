

% trange - defines the range for the block average [tPre tPost dt]. If dt
%           defined, time series are interpolated prior to 

function [data_yavg, data_yavgstd] = fUS_GLM_scratch(data, t, stim, Aaux, tIncAuto, trange, glmSolveMethod, idxBasis, paramsBasis, driftOrder)

    %%%% Get characteristics of the timing
    dt = t(2) - t(1); % Time step
    fq = 1/dt;        % Sampling/measurement frequency
    nPre = round(trange(1)/dt); % Which index for the block average Pre point
    nPost = round(trange(2)/dt);% Which index for the block average Post point
    % nTpts = size(y,1); % # of time points
    nTpts = length(t); % # of time points
    tHRF = (1*nPre*dt:dt:nPost*dt)';
    ntHRF = length(tHRF); % # of points in the HRF
    nT = length(t); % should be the same as nTpts (?)

    % % % % % % % %%%% Update stimulus vectors
    % % % % % % % % Note: it's assumed that the stimulus vector has the same length and
    % % % % % % % % timesteps of the data vector.
    % % % % % % % % stim represents the onsets of the active stimulus periods.
    % % % % % % % stimStates = stim;
    % % % % % % % stimStates(stimStates > 0) = stimStates(stimStates > 0) ./ stimStates(stimStates > 0);
    % % % % % % % stimStates = logical(stimStates); % stimStates: a logical version of the stim vector. True = stim on, false = stim off
    % % % % % % % stimAmps = stim; % stimAmps: a quantitative version of the stim vector. The amplitude of the stimulus is preserved.
    % % % % % % % 
    % % % % % % % %%%%%%%%%%%%%%%%
    % % % % % % % % Prune good stim, generate onset matrix
    % % % % % % % %%%%%%%%%%%%%%%%
    % % % % % % % % Get only indices of conditions with any stimStates that are 1
    % % % % % % % % lstCond = find(sum(stimStates == 1, 1) > 0);
    % % % % % % % lstCond = find(stimStates == true); % Indices at the active stim onsets
    % % % % % % % % nCond = length(lstCond); % # of stim onsets
    % % % % % % % % nTrials = zeros(nCond, 1); % # of trials??......
    % % % % % % % nCond = 1; % Only one condition
    % % % % % % % nTrials = length(lstCond);
    % % % % % % % 
    % % % % % % % onset = zeros(nT, nCond);
    % % % % % % % avg_pulses = {};
    % % % % % % % for iCond = 1:nCond
    % % % % % % %     % lstT = find(stimStates(:, lstCond(iCond)) == 1);  % Indices of stims enabled (== 1)
    % % % % % % %     % lstp = find((lstT+nPre) >= 1 & (lstT+nPost) <= nTpts);  % Indices of stims not clipped by signal
    % % % % % % %     lstT = find(stimStates == true);  % Indices of stims enabled (== 1)
    % % % % % % %     lstp = find((lstT+nPre) >= 1 & (lstT+nPost) <= nTpts);  % Indices of stims not clipped by signal
    % % % % % % %     lst = lstT(lstp); % Final list of stim onset indices to use
    % % % % % % %     % nTrials(iCond) = length(lst);
    % % % % % % % 
    % % % % % % %     % Generate basis boxcars of stim amplitude and duration
    % % % % % % %     starts = lst+nPre; % Get the start indices of each trial's "stim" --> include the Pre offset
    % % % % % % %     % if ~isempty(stim(lstCond(iCond)))
    % % % % % % %         durations = stim(lstCond(iCond)).data(:, 2);
    % % % % % % %         amplitudes = stim(lstCond(iCond)).data(:, 3);
    % % % % % % %         avg_pulses{iCond} = ones(round(mean(durations) / dt), 1); %#ok<AGROW>
    % % % % % % %         for i = 1:length(starts)
    % % % % % % %             if idxBasis == 1  % Gaussian has no duration T (yet)
    % % % % % % %                pulse_duration = 1; 
    % % % % % % %             else
    % % % % % % %                pulse_duration = round(durations(i) / dt); 
    % % % % % % %             end
    % % % % % % %             pulse = (amplitudes(i) / pulse_duration) * ones(pulse_duration, 1);
    % % % % % % %             onset(starts(i):starts(i) + pulse_duration - 1, iCond) = onset(starts(i):starts(i) + pulse_duration - 1, iCond) + pulse;
    % % % % % % %         end
    % % % % % % %     % end
    % % % % % % % end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct the basis functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % nB: # of basis functions to use
    % tbasis = "temporal" basis, probably
    switch idxBasis
        case 1
            % Gaussians
            gms = paramsBasis(1);  % Mean
            gstd = paramsBasis(2); % Standard deviation
            
            nB = floor((trange(2)-trange(1)) / gms) - 1;
            tbasis = zeros(ntHRF, nB); % temp (?) basis matrix: # of HRF time points by # of basis functions
            for b=1:nB % go through each basis function
                tbasis(:, b) = exp(-(tHRF-(trange(1)+b*gms)).^2/(2*gstd.^2));
                tbasis(:, b) = tbasis(:,b)./max(tbasis(:,b));
            end
        
        case 2
            % Modified Gamma
            % if length(paramsBasis)==2
                nConc = 1; % Should only need 1 set of basis functions. Theirs should stand for multiple chromophores (HbO, HbR)
            % else
            %     nConc = 2;
            % end
            
            nB = 1; % 1 basis
            tbasis = zeros(ntHRF, nB, nConc);
            for iConc = 1:nConc
                tau = paramsBasis((iConc-1)*2+1);
                sigma = paramsBasis((iConc-1)*2+2);
                
                tbasis(:,1,iConc) = (exp(1)*(tHRF-tau).^2/sigma^2) .* exp( -(tHRF-tau).^2/sigma^2 );
                lstNeg = find(tHRF<0); % List of indices where the HRF time is negative 
                tbasis(lstNeg, 1, iConc) = 0; % ^ Set those points' basis values to 0
                
                if tHRF(1)<tau
                    tbasis(1:round((tau-tHRF(1))/dt),1,iConc) = 0; % Set the basis values to 0 when ...
                end
                
            end
        
        case 3
            % Modified Gamma and Derivative
            % if length(paramsBasis)==2
                nConc = 1;
            % else
            %     nConc = 2;
            % end
            
            nB = 2;
            tbasis = zeros(ntHRF, nB, nConc);
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
            % if length(paramsBasis)==2
                nConc = 1;
            % else
            %     nConc = 2;
            % end
            
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
    % dA=zeros(nT, nB*nCond, 2);
    dA=zeros(nT, nB*nCond);
    % for iConc = 1:2
    for iConc = 1
        iC = 0;
        for iCond=1:nCond % Go through each condition (should be just 1 for us)
            for b=1:nB % Go through each basis
                iC = iC + 1;

                % Convolve the basis functions with the boxcars to get the
                % actual bases for GLM
                if size(tbasis,3)==1
                    clmn = conv(onset(:,iCond),tbasis(:,b));
                else
                    clmn = conv(onset(:,iCond),tbasis(:,b,iConc));
                end
                clmn = clmn(1:nT);
                dA(:, iC, iConc) = clmn;
                beta_label{b + (iCond-1)*nB} = ['Cond' num2str(iCond)];
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Expand design matrix nth order polynomial for drift correction
    % rescale polynomial to avoid bad conditionning
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xDrift = ones(nT, driftOrder);
    for ii = 2:(driftOrder+1) % ii = n + 1...
        xDrift(:, ii) = ([1:nT]').^(ii-1); % Create each polynomial
        xDrift(:, ii) = xDrift(:, ii) / xDrift(end, ii);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Expand design matrix with Aaux
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nAux = size(Aaux, 2);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Expand design matrix for Motion Correction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     if flagMotionCorrect==1
    %         idxMA = find(diff(tInc)==1);  % number of motion artifacts
    %         if isempty(idxMA)
    nMC = 0;
    Amotion = [];
    %         else
    %             nMA = length(idxMA);
    %             nMC = nMA+1;
    %             Amotion = zeros(nT,nMC);
    %             Amotion(1:idxMA(1),1) = 1;
    %             for ii=2:nMA
    %                 Amotion((idxMA(ii-1)+1):idxMA(ii),ii) = 1;
    %             end
    %             Amotion((idxMA(nMA)+1):end,end) = 1;
    %         end
    %     else
    %         nMC = 0;
    %         Amotion = [];
    %     end
    lstInc = find(tInc==1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Final design matrix (A)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dummy = size(beta_label,2);

    for iConc = 1:nConc
        A(:, :, iConc)=[dA(:, :, iConc) xDrift Amotion];
    end

    nCh = 1; % # of channels (1 for ultrasound)

    % Exit if not enough data to analyze; the 3 here is arbitrary.
    % Certainly needs to be larger than 1
    if length(lstInc) < 3*size(A, 2) || nCond==0
        warning('Not enough data to find a solution')
        yavg    = zeros(ntHRF, nCh, 3,nCond);
        yavgstd = zeros(ntHRF, nCh, 3,nCond);
        ysum2   = zeros(ntHRF, nCh, 3,nCond);
        return
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SOLVE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tb = zeros(nB, nCh, nCond);
    % b = zeros(driftOrder+1+nAux,nCh,2);
    yavg    = zeros(ntHRF, nCh, nCond);
    yavgstd = zeros(ntHRF, nCh, nCond);
    ysum2   = zeros(ntHRF, nCh, nCond);
    yresid  = zeros(nT, nCh);
    ynew    = zeros(nT, nCh);
    yR      = zeros(nCh);
    foo     = zeros(nB*nCond+driftOrder+1+nAux+nMC,nCh,2); % 4 extra for 3rd order drift + nAux


    % check if the matrix is well conditionned
    ATA=At(lstInc,:)'*At(lstInc,:); % what is lstInc doing???????
    rco = rcond(full(ATA)); % Reciprocal condition
    if rco<10^-14 && rco>10^-25
        fprintf('Design matrix is poorly scaled...(RCond=%e)\n', rco);
    elseif rco<10^-25
        fprintf('Design matrix is VERY poorly scaled...(RCond=%e), cannot perform computation\n', rco);
        yavg = permute(yavg,[1 3 2 4]);
        yavgstd = permute(yavgstd,[1 3 2 4]);
        ysum2 = permute(ysum2,[1 3 2 4]);
        ynew = y;
        yresid = zeros(size(y));
        
        foo = nTrials;
        nTrials = zeros(1,size(stimAmps,2));
        nTrials(lstCond) = foo;
        
        foo = yavg;
        yavg = zeros(size(foo,1),size(foo,2),size(foo,3),size(stimAmps,2));
        yavg(:,:,:,lstCond) = foo;
        
        foo = yavgstd;
        yavgstd = zeros(size(foo,1),size(foo,2),size(foo,3),size(stimAmps,2));
        yavgstd(:,:,:,lstCond) = foo;
        
        foo = ysum2;
        ysum2 = zeros(size(foo,1),size(foo,2),size(foo,3),size(stimAmps,2));
        ysum2(:,:,:,lstCond) = foo;
        
        beta = [];
        return
    end
    
    % Compute pseudo-inverse and deconvolve
    if glmSolveMethod==1 % ~flagUseTed
        pinvA=ATA\At(lstInc,:)';
        foo = [];
        ytmp = y(lstInc,conc,lstML);
        foo(:,lstML,conc)=pinvA*squeeze(ytmp);
    elseif glmSolveMethod==2
        % Use the iWLS code from Barker et al
        foo = [];
        ytmp = y(lstInc,conc,lstML);
        for chanIdx=1:length(lstML)
            ytmp2 = y(lstInc,conc,lstML(chanIdx));
            [dmoco, beta, tstat(:,lstML(chanIdx),conc), pval(:,lstML(chanIdx),conc), sigma, CovB(:,:,lstML(chanIdx),conc), dfe, w, P, f] = ar_glm_final(squeeze(ytmp2),At(lstInc,:), round(fq*2));

            foo(:,lstML(chanIdx),conc)=beta;
            ytmp(:,1,chanIdx) = dmoco; %We also need to keep my version of "Yvar" and "Bvar"                    
            
            yvar(:,lstML(chanIdx),conc)=sigma.^2;
            bvar(:,lstML(chanIdx),conc)=diag(CovB(:,:,lstML(chanIdx),conc));  %Note-  I am only keeping the diag terms.  This lets you test if beta != 0,
            %but in the future the HOMER-2 code needs to be modified to keep the entire cov-beta matrix which you need to test between conditions e.g. if beta(1) ~= beta(2)
            % reply to above comment from DAB: Now directly using CovB for the contrast calculations. MAY
        end
    end
    
    % Solution
    for iCond=1:nCond
        tb(:,lstML,conc,iCond)=foo([1:nB]+(iCond-1)*nB,lstML,conc);
        %                yavg(:,lstML,conc,lstCond(iCond))=tbasis*tb(:,lstML,conc,lstCond(iCond));
        if size(tbasis,3)==1
            yavg(:,lstML,conc,iCond)=tbasis*tb(:,lstML,conc,iCond);
        else
            yavg(:,lstML,conc,iCond)=tbasis(:,:,conc)*tb(:,lstML,conc,iCond);
        end
        if idxBasis > 1
            for iML = transpose(lstML)
                convolved = conv(yavg(:, iML, conc, iCond), avg_pulses{iCond});
                yavg(:, iML, conc, iCond) = convolved(1:size(yavg, 1));  % Truncate convolution
            end
        end
    end
    
    % Reconstruct y and yresid (y is obtained just from the HRF) and R
    yresid(lstInc,conc,lstML) = ytmp - permute(At(lstInc,:)*foo(:,lstML,conc),[1 3 2]);
    ynew(lstInc,conc,lstML) = permute(dA(lstInc,:,conc)*foo(1:(nB*nCond),lstML,conc),[1 3 2]) + yresid(lstInc,conc,lstML);
    
    yfit = permute(At(lstInc,:)*foo(:,lstML,conc),[1 3 2]);
    for iML=1:length(lstML)
        yRtmp = corrcoef(ytmp(:,1,iML),yfit(:,1,iML));
        yR(lstML(iML),conc) = yRtmp(1,2);
    end
    
    % Get error
    if glmSolveMethod==1 %  OLS  ~flagUseTed
        pAinvAinvD = diag(pinvA*pinvA');
        yest(:,lstML,conc) = At * foo(:,lstML,conc);
        yvar(1,lstML,conc) = sum((squeeze(y(:,conc,lstML))-yest(:,lstML,conc)).^2)./(size(y,1)-1); % check this against eq(53) in Ye2009
        for iCh = 1:length(lstML)
            
            % GLM stats for each condition
            bvar(:,lstML(iCh),conc) = yvar(1,lstML(iCh),conc) * pAinvAinvD;
            tval(:,lstML(iCh),conc) =  foo(:,lstML(iCh),conc)./sqrt(bvar(:,lstML(iCh),conc));
            pval(:,lstML(iCh),conc) = 1-tcdf(abs(tval(:,lstML(iCh),conc)),(size(y,1)-1));
            %
            
            % GLM stats for contrast between conditions, given a c_vector exists
            if nCond > 1
                if (sum(abs(c_vector)) ~= 0) && (size(c_vector,2) == nCond)
                    
                    if ~exist('cv_extended','var') == 1
                        cv_dummy = [];
                        for m = 1:nCond
                            cv_dummy = [cv_dummy ones(1,nB)*c_vector(m)];
                        end
                        cv_extended = [cv_dummy zeros(1,size(beta_label,2)-size(cv_dummy,2))];
                    end
                    
                    tval_contrast(:,lstML(iCh),conc) = cv_extended * foo(:,lstML(iCh),conc)./sqrt(cv_extended * (pinvA*pinvA') * yvar(:,lstML(iCh),conc) * cv_extended');
                    pval_contrast(:,lstML(iCh),conc) = 1-tcdf(abs(tval_contrast(:,lstML(iCh),conc)),(size(y,1)-1));
                end
            end
            %
            
            
            for iCond=1:nCond
                if size(tbasis,3)==1
                    yavgstd(:,lstML(iCh),conc,iCond) = diag(tbasis*diag(bvar([1:nB]+(iCond-1)*nB,lstML(iCh),conc))*tbasis').^0.5;
                else
                    yavgstd(:,lstML(iCh),conc,iCond) = diag(tbasis(:,:,conc)*diag(bvar([1:nB]+(iCond-1)*nB,lstML(iCh),conc))*tbasis(:,:,conc)').^0.5;
                end
                ysum2(:,lstML(iCh),conc,iCond) = yavgstd(:,lstML(iCh),conc,iCond).^2 + nTrials{iBlk}(iCond)*yavg(:,lstML(iCh),conc,iCond).^2;
            end
        end
        
    elseif glmSolveMethod==2  % iWLS
        
        yest(:,lstML,conc) = At * foo(:,lstML,conc);
        for iCh = 1:length(lstML)
            
            % GLM stats for contrast between conditions, given a c_vector exists
            if nCond > 1
                if (sum(abs(c_vector)) ~= 0) && (size(c_vector,2) == nCond)
                    
                    if ~exist('cv_extended','var') == 1
                        cv_dummy = [];
                        for m = 1:nCond
                            cv_dummy = [cv_dummy ones(1,nB)*c_vector(m)];
                        end
                        cv_extended = [cv_dummy zeros(1,size(At,2)-size(cv_dummy,2))];
                    end
                    tval_contrast(:,lstML(iCh),conc) = cv_extended * foo(:,lstML(iCh),conc)./sqrt(cv_extended * squeeze(CovB(:,:,lstML(chanIdx),conc))* cv_extended');
                    pval_contrast(:,lstML(iCh),conc) = 1-tcdf(abs(tval_contrast(:,lstML(iCh),conc)),(size(y,1)-1));
                end
            end
            %
            
            for iCond=1:nCond
                if size(tbasis,3)==1
                    yavgstd(:,lstML(iCh),conc,iCond) = diag(tbasis*diag(bvar([1:nB]+(iCond-1)*nB,lstML(iCh),conc))*tbasis').^0.5;
                else
                    yavgstd(:,lstML(iCh),conc,iCond) = diag(tbasis(:,:,conc)*diag(bvar([1:nB]+(iCond-1)*nB,lstML(iCh),conc))*tbasis(:,:,conc)').^0.5;
                end
                ysum2(:,lstML(iCh),conc,iCond) = yavgstd(:,lstML(iCh),conc,iCond).^2 + nTrials{iBlk}(iCond)*yavg(:,lstML(iCh),conc,iCond).^2;
            end
        end
        
    end
end