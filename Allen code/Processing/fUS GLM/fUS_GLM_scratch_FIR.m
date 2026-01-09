
% ti: (timing data) a struct with the fields:
%       - t: [# time points x 1] vector of timestamps [probably in seconds]
%       - stimAmps: [# time points x 1] vector of stim amplitudes (0 = off, > 0 = on)
%       - tOnsets: [# stims x 1] vector of time points corresponding to stim starts/onsets
%       - (not implemented) stimOnsets: [# time points x 1] vector of stim starts/onsets (0 = off, 1 = start)
% trange - defines the range for the block average [tPre tPost dt]. If dt
%           defined, time series are interpolated prior to 

function [data_yavg, data_yavgstd] = fUS_GLM_scratch_FIR(data, ti, trange, FIR_order, glmSolveMethod, paramsBasis)

    t = ti.t; % Get the time vector from the input struct

    %%%% Get characteristics of the timing
    % dt = t(2) - t(1); % Time step
    dt = median(diff(t)); % Time step; use mean of all the timesteps in case the sampling is not uniform
    fq = 1/dt;        % Sampling/measurement frequency
    nPre = round(trange(1)/dt);  % Which index for the block average Pre point
    nPost = round(trange(2)/dt); % Which index for the block average Post point
    % nTpts = size(y,1); % # of time points
    nTpts = length(t); % # of time points
    tHRF = (nPre * dt:dt:nPost * dt)'; % Time points for the HRF
    ntHRF = length(tHRF); % # of points in the HRF
    nT = length(t); % should be the same as nTpts (?)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % My version of the stim info
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct the basis functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tbasis = ti.stimAmps; % is this right??? does it have the correct time stamps?

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct design matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A = repmat(tbasis, 1, FIR_order); % ???

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Final design matrix (A)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nRegs = size(beta_label, 2); % A dummy variable that keeps track of the number of regressors used (previously called 'dummy')

    
    % Update beta (weight) labels
    for ixDrift = 1:size(xDrift, 2) % Label each drift (order)
        beta_label{ixDrift + nRegs} = "xDrift: order " + num2str(ixDrift - 1);
    end
    nRegs = size(beta_label, 2);
    for iAaux = 1:size(Aaux, 2)    % Label each auxiliary regressor
        beta_label{iAaux + nRegs} = "Aux # " + num2str(iAaux);
    end
    nRegs = size(beta_label, 2);
    for iAmotion = 1:size(Amotion, 2)
        beta_label{iAmotion + nRegs} = "Motion # " + num2str(iAmotion);
    end
    nRegs = size(beta_label, 2);

    nCh = 1; % # of channels (1 for ultrasound)

    % Exit if not enough data to analyze; the 3 here is arbitrary.
    % Certainly needs to be larger than 1
    if length(lstInc) < 3*size(A, 2) || nCond==0 % If there are much fewer included time points than the number of basis functions...
        warning('Not enough data to find a solution')
        yavg    = zeros(ntHRF, nCh, nCond);
        yavgstd = zeros(ntHRF, nCh, nCond);
        ysum2   = zeros(ntHRF, nCh, nCond);
        return
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SOLVE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tb = zeros(nB, nCh, nCond);
    yavg    = zeros(ntHRF, nCh, nCond);
    yavgstd = zeros(ntHRF, nCh, nCond);
    ysum2   = zeros(ntHRF, nCh, nCond);
    yresid  = zeros(nT, nCh);
    ynew    = zeros(nT, nCh);
    yR      = zeros(nCh);
    foo     = zeros(nB*nCond + (driftOrder + 1) + 1 + nAux + nMC, nCh); % original comment: 4 extra for 3rd order drift + nAux


    % check if the matrix is well conditionned
    ATA = A(lstInc, :)' * A(lstInc, :); % lstInc keeps the "non-motion" points (if that's inputted)
    rco = rcond(full(ATA)); % Reciprocal condition
    if rco < 10^-14 && rco > 10^-25
        fprintf('Design matrix is poorly scaled...(RCond=%e)\n', rco);
    elseif rco < 10^-25
        fprintf('Design matrix is VERY poorly scaled...(RCond=%e), cannot perform computation\n', rco);
        yavg = permute(yavg,[1 3 2 4]);
        yavgstd = permute(yavgstd,[1 3 2 4]);
        ysum2 = permute(ysum2,[1 3 2 4]);
        ynew = y;
        yresid = zeros(size(y));
        
        % ????
        foo = nTrials;
        nTrials = zeros(1, size(ti.stimAmps, 2));
        nTrials(lstCond) = foo;
        
        foo = yavg;
        yavg = zeros(size(foo, 1), size(foo, 2), size(foo, 3), size(ti.stimAmps, 2));
        yavg(:, :, :,lstCond) = foo;
        
        foo = yavgstd;
        yavgstd = zeros(size(foo, 1), size(foo, 2), size(foo, 3), size(ti.stimAmps, 2));
        yavgstd(:, :, :, lstCond) = foo;
        
        foo = ysum2;
        ysum2 = zeros(size(foo, 1),size(foo, 2), size(foo, 3), size(ti.stimAmps, 2));
        ysum2(:, :, :, lstCond) = foo;
        
        beta = [];
        return
    end
    
    % Compute pseudo-inverse and deconvolve
    if glmSolveMethod == 1 % Least squares (original comment: ~flagUseTed)

        pinvA = ATA \ A(lstInc, :)'; % Solve for weights
        foo = [];
        ytmp = y(lstInc, conc, lstML); % ??????????????????????????????????????????????????????????
        foo(:,lstML,conc) = pinvA*squeeze(ytmp);
    elseif glmSolveMethod == 2 % Iterative weighted least squares
        % Use the iWLS code from Barker et al
        foo = [];
        ytmp = y(lstInc,conc,lstML); %..............................
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