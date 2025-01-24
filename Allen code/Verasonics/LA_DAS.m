% Linear Array Delay and Sum
% IQ: DAS beamformed data
% pixSpacing: [z pixel spacing (m), x pixel spacing (m)]
function [IQ, pixSpacing] = LA_DAS(RcvData, P, zPixSpacing)
    %% changes
    % 9/12/24: now properly accounts for a nonzero startDepth
    % 11/7/24: now a function
    %% potential to do/issues
    
    % I feel like I'm still missing something that improves quality
    
    % add band pass filter?
    %%
    useGain = 1;
    % Currently it's set to use the element sensitivity as a weight
        
    %% Start parpool
    % https://www.mathworks.com/matlabcentral/answers/91744-how-can-i-check-if-matlabpool-is-running-when-using-parallel-computing-toolbox
    
    pp = gcp('nocreate');
    if isempty(pp)
        % There is no parallel pool
        parpool LocalProfile1
    
    end
    
     %% Load variables into the function space so you don't waste more time accessing the structure each time
    % have to initialize like this or it doesn't work...
    tw_peak = [];
    angles = [];
    endDepth = [];
    endDepthMM = [];
    Event = [];
    fps_target = [];
    % flow_v_mm_s = [];
    maxAcqLength_adjusted = [];
    maxAngle = [];
    Media = [];
    na = [];
    nf = [];
    numElements = [];
    Receive = [];
    Resource = [];
    SeqControl = [];
    startDepth = [];
    startDepthMM = [];
    TGC = [];
    Trans = [];
    TW = [];
    TX = [];
    wl = [];
    L = [];
    numFramesPerBuffer = [];
    
    assignStructVars(P)
    
    nf = numSubFrames;
%     if exist(P.numSubFrames)
%         nf = numSubFrames;
%     elseif exist(P.numFramesPerBuffer)
%         nf = numFramesPerBuffer;
%     else
%         error('something wrong with the P struct, missing num frames variable')
%     end
    tw_peak = TW.peak;
    
    %%
    tstart = clock;
    
    if iscell(RcvData)
        r = RcvData{1, 1}; % use this for the simulation where data is saved together at the end
    else
        r = RcvData;         % use this for continuous acquisition (ctsacq) when data is saved after each set of frames
    end
    
    %     r_old = r;
    % ff = 200;
    % r = r(:, :, 1:ff);
    % nf = ff;
    
    r = r(:, Trans.Connector, :); % Reorganize data according to array element-to-channel map
    
    s = size(r);

    % temp fix 1/9/25
    s(1) = Receive(1).endSample*na;
    
    A = Media.attenuation; % dB/cm/MHz
    A_per_m = A * Trans.frequency * 100; % dB/m
    
    c = Resource.Parameters.speedOfSound;
    tLensCorrection = Trans.lensCorrection;
    %     t_correction = 0; % Transducer lens correction
    t_correction = 2 * (tLensCorrection * wl) / c; % Trans.lensCorrection is the one way delay in wavelengths through lens
    d_correction_1way = (tLensCorrection * wl);     % use further samples to account for delay.......
    
    element_sens_cutoff = 0.6; %%%%%%%%%%%%%%%%%%%%%%%%%
    %     samplesPerWL = 4;
    samplesPerWL = Receive(1).samplesPerWave; % Might not be exactly 4 due to their sampling
    
    timePerSample = 1/(Trans.frequency * 1e6)/samplesPerWL;
    distPerSample = wl/samplesPerWL;
    
    %     if length(s) == 3 % 1D array/2D data
    
    %%%%%%%%%%%%%%%%%%%%
    % reorganize RcvData into (# z samples, # x samples, # angles, # frames)
    nzs_nopad = s(1)/na;
    nxs_nopad = s(2);
    d = zeros(nzs_nopad, nxs_nopad, na, nf);
    % d = zeros(Receive(1).endSample, s(2), na, nf);
    parfor f = 1:nf
        for a = 1:na
            d(:, :, a, f) = r(Receive(a).startSample:Receive(a).endSample, :, f);
        end
    end
    % d_old = d; %%%%% for testing purposes

    d = hilbert(d);
    
    d_pad = zeros(round(startDepth * samplesPerWL), nxs_nopad, na, nf);
    d_pad_size = size(d_pad);
    d = [d_pad; d; d_pad];
    d = squeeze(d);
    
    % d = d(:, :, 2, :); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    d_size = size(d);
    
    % imagesc(squeeze(d(:, :, 11, 1)))
    
    
    nzs = d_size(1);
    nxs = d_size(2);
    
    clear r % save space
    
    %%
    % pixsizez = wl/2;
    pixsizez = zPixSpacing;
%     pixsizex = wl;
    pixsizex = Trans.spacingMm/1e3;
    zsize = (endDepth - startDepth) * wl;    % region's z span in m
    xsize = numElements.*Trans.spacing * wl; % region's x span in m
            znumpix = ceil(zsize / pixsizez);         % # rows
    %         znumpix = 712;
    %         xnumpix = ceil(xsize / pixsizex);         % # cols
    
%     znumpix = nzs_nopad; % Full z sample size
    xnumpix = nxs_nopad;
    
    % znumpix = nzs; % Full z sample size
    % xnumpix = nxs;
    pixSpacing = [pixsizez, pixsizex];
    %% plot raw waveform from one element
    % 
    % figure
    % plot(d(:, round(Media.MP(1, 1) + numElements/2), 1, 1))
    % hold on
    % plot(d(:, round(Media.MP(1, 1) + numElements/2), 3, 1))
    % plot(d(:, round(Media.MP(1, 1) + numElements/2), 5, 1))
    % hold off
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    arrayLength = numElements * Trans.spacingMm / 1e3; % transducer length in m
    elementx = (wl .* Trans.ElementPos(:, 1))';
    x_end = max(elementx);
    
    d_sum = zeros(znumpix, xnumpix, na, nf);
    
    d_delayed_temp = zeros(numElements, nf); % get the value at a single sample point, supposed to be the max of the envelope
                                       % (# x samples, # frames)
    
    % make this more efficient?
    frame_ind_temp = zeros(numElements * nf * na, 1);
    for f = 1:nf
        frame_ind_temp((f-1)*numElements * na + 1 : f * na * numElements) = f;
    end
    
    angles_ind_temp = zeros(numElements * na, 1);
    for a = 1:na
        angles_ind_temp((a-1)*numElements + 1 : a * numElements) = a;
    end
    angles_ind_temp = repmat(angles_ind_temp, nf, 1);
    
    %%
    parfor zi = 1:znumpix
    % for zi = 877
    % for zi = 1:znumpix
    % for zi = znumpix
    % for zi = round(zl/(endDepth) * znumpix)
    % for zi = round(1.038688524590164e+02)
    %     z_level_pix = zi / znumpix * zsize;      % z value of zi in m
    
    %     z_level_pix = startDepth * wl + zi / znumpix * zsize;  % z value of zi in m. Start z recon at startDepth * wl
    %     z_level_pix = zi / znumpix * zsize;  % z value of zi in m.
        % I think the below one is more correct, but the above one works better
            z_level_pix = startDepth * wl + (zi - 1/2) / znumpix * zsize;  % z value of zi in m. Start z recon at startDepth * wl
    
        for xi = 1:xnumpix
    %     for xi = 64
    %     for xi = 1
    %     for xi = xnumpix
    %     for xi = round(xl)
    %     for xi = round(Media.MP(1, 1) + numElements/2)
    %     for xi = 25
    
    %         x_level_pix = (xi - xnumpix/2) / xnumpix * xsize;  % x value of xi in m. x = 0 is center of transducer array
            % I think the below one is more correct, but the above one works
            % better (nvm maybe????)
    %         x_level_pix = (xi - xnumpix/2 - 1/2) / xnumpix * xsize;  % x value of xi in m. x = 0 is center of transducer array
            
            % change 9/24/24
            x_level_pix = (xi - (xnumpix + 1)/2) / xnumpix * xsize;  % x value of xi in m. x = 0 is center of transducer array
    
    %         traveldistTX_temp = z_level_pix .* cos(angles) + (sign(angles) .* x_end - x_level_pix) .* sin(angles);    % travel dist from TX to pixel depends on z and angle
    %         traveldistTX_temp = (z_level_pix - startDepth*wl) .* cos(angles) + (sign(angles) .* x_end - x_level_pix) .* sin(angles);    % travel dist from TX to pixel depends on z and angle
            
            %%%%
            
    %         xa = x_level_pix + (z_level_pix - startDepth * wl) .* tan(angles);
            
    %         za = startDepth * wl;
    %         xa = elementx + za .* tan(angle_i);
            lcTX = d_correction_1way;
    %         if startDepth > tLensCorrection
    %             lcTX = 0;
    %         else
    %             lcTX = d_correction_1way;
    %         end
    
            % change 9/24/24 - disable the fliplr on ind
            traveldistTX_temp = z_level_pix .* cos(angles) + (x_level_pix + sign(angles) .* x_end) .* sin(angles) + 1 * d_correction_1way;
    
            traveldistTX_temp = repmat(traveldistTX_temp, numElements, 1);
    
            % pixel to element receive distance
            traveldistRX_temp = (sqrt(z_level_pix ^ 2 + (elementx - x_level_pix).^2) + d_correction_1way)'; % RX travel dist from pixel to each transducer element i
    
    %         temp_td = traveldistTX_temp + repmat(traveldistRX_temp, 1, na) + wl/2;
            % change 9/24/24
            temp_td = traveldistTX_temp + repmat(traveldistRX_temp, 1, na) + wl * (tw_peak);
    
    %                 temp_td = traveldistTX_temp + repmat(traveldistRX_temp, 1, na) + wl * ((TW.peak) + 0.25);
    
    %           temp_td = traveldistTX_temp + repmat(traveldistRX_temp, 1, na);
    
    %         temp_td = traveldistTX_temp + repmat(traveldistRX_temp, 1, na) + 2 * d_correction_1way;
    %         disp(max(temp_td))
    %         temp_td = traveldistTX_temp + repmat(traveldistRX_temp, 1, na) + 2 * d_correction_1way - startDepth*wl;
    
        
            % do I want to unaccount for the angle if the pixel is out
            % of the beam
    
            %%%%%%%%%%%%%%% interp???????
            % round seems to work better than floor
    
            ind = round(temp_td ./ distPerSample - startDepth*samplesPerWL);
    
            %%%%%%%%%%%%%%%%%%%%%%
            ind(ind > nzs) = nzs;
            %%%%%%%%%%%%%%%%%%%%%%
            
            ia = repmat(reshape(ind, nxs*na, 1), nf, 1);
            ib = repmat((1:numElements)', nf*na, 1);
            ind_linear = sub2ind(d_size, ia, ib, angles_ind_temp, frame_ind_temp);
    
            d_delayed_temp = d(ind_linear); % vector that stacks the delayed d value from each frame
            
            dd = zeros(nxs, na, nf); 
            for f = 1:nf
                for a = 1:na
                    dd(:, a, f) = d_delayed_temp((f-1)*na * numElements + (a-1)*numElements + 1 : (f-1)*na*numElements + (a)*numElements);
    
                end
            end
            dd = squeeze(dd);
    
            % Account for element angle sensitivity 
             % angle from a pixel to each element
            angle_i = atan((x_level_pix - elementx) ./ z_level_pix); % I think atan and not atan2 is right here. ccw is +
    
            % interpolate provided element vs. angle sensitivity curve
    %         element_sens_i = interp1(linspace(-pi/2, pi/2, length(Trans.ElementSens)), Trans.ElementSens, angle_i);
    %         element_sens_mask = element_sens_i > element_sens_cutoff;
    % 
    %         dd(repmat(~element_sens_mask', 1, na, nf)) = 0;
    
            % change 10/21/24: weighting
            element_sens_i = interp1(linspace(-pi/2, pi/2, length(Trans.ElementSens)), Trans.ElementSens, angle_i)';
            dd = dd .* repmat(element_sens_i, 1, na, nf);
     
            % change 10/21/24
            if useGain
                dd = dd .* 10.^((temp_td .* -A_per_m/1)  ./ 20);
            end
    
            d_sum(zi, xi, :, :) = squeeze(sum(dd, 1)); % sum delayed signal values, across elements
    
        end
        
    end
    
    
    %%
        tend = clock;
        disp(etime(tend, tstart))
    
        % Calculate IQ
        IQ = d_sum;
%         IQ = hilbert(d_sum); % operates on the columns (z) of the d_sum matrix
    %     env = abs(IQ);       % Absolute value of the Hilbert transform is the envelope of the original signal
    
    
        % imagesc(squeeze(real(IQ(:, :, 1, 1))))
    
    %         % Single angle image
    %         figure(1)
    %         imagesc(squeeze(env(:, :, 1, 1)))
    
    %                 imagesc(squeeze(d_sum(:, :, 1, 1)))
    
        % Get phase of each pixel in the image
    %     phase = atan2(imag(IQ), real(IQ));     % Need atan2 for full 0 - 2pi range!!
    
    %         plot(phase(:, 64, 1))
    % plot(d_delayed{560, 64}(:, 1, 1))
    
    
    
    % %% Coherently sum across angles and plot
    % 
    % % % ap = env + 1i.*phase;
    % % % ap_sum = sum(ap, 3);
    % % % figure(3)
    % % % ap_sum_intensity = squeeze(abs(ap_sum));
    % % % imagesc(ap_sum_intensity(:, :, 1))
    % % % title('Probably wrong')
    % 
    % apt = squeeze(sum(IQ, 3));
    % figure
    % apt_sum_intensity = abs(apt);
    % imagesc(apt_sum_intensity(:, :, 1))
    % title('Coherently summed image, frame 1')
    % % ax=gca;
    % % yticks./znumpix.*endDepth.*wl
    % % 
    % % 
    % % Plot Verasonics' reconstructed image for one angle
    % % figure(2)
    % % iqv = ImgData{1, 1};
    % 
    % % imagesc(iqv(:, :, 1))
    % 
    % % Verasonics coherent sum
    % IQv = squeeze(IData{1, 1} + 1i.*QData{1, 1});
    % figure
    % IQv_sum = squeeze(sum(IQv, 3));
    % [znumpixv, xnumpixv, ~] = size(IQv_sum);
    % abs_IQv_sum = abs(IQv_sum);
    % imagesc(abs(abs_IQv_sum(:, :, 1)))
    % title('Verasonics sum')
    % % 
    %% 1D PSFs for our reconstruction after coherent summation
    % figure(6)
    % % plot(apt_sum_intensity(round((Media.MP(4, 3) - startDepth) * samplesPerWL), :, 1))
    % plot(apt_sum_intensity(697, :, 1), 'LineWidth', 1.5)
    % title(strcat('Lateral PSF at z = ', num2str(Media.MP(4, 3)), ' wavelengths'))
    % xlim([1, xnumpix])
    % xlabel('x pixel')
    % ylabel('Intensity')
    % % 
    % figure(7)
    % plot(apt_sum_intensity(:, 64, 1), 'LineWidth', 1.5)
    % title(strcat('Axial PSF at x = ', num2str(64 - numElements/2), ' wavelengths'))
    % % title('Axial PSF')
    % % xlim([d_pad_size(1), znumpix])
    % xlim([1, znumpix])
    % xlabel('z pixel')
    % ylabel('Intensity')
    % 
    % figure(8)
    % plot(apt_sum_intensity(:, 104, 1), 'LineWidth', 1.5)
    % title(strcat('Axial PSF at x = ', num2str(104 - numElements/2), ' wavelengths'))
    % % title('Axial PSF')
    % % xlim([d_pad_size(1), znumpix])
    % xlim([1, znumpix])
    % xlabel('z pixel')
    % ylabel('Intensity')
    
    %% 1D PSFs for Verasonics reconstruction after coherent summation
    
    % figure
    % % plot(apt_sum_intensity(round((Media.MP(4, 3) - startDepth) * samplesPerWL), :, 1))
    % plot(abs_IQv_sum(119, :, 1), 'LineWidth', 1.5)
    % title(strcat('Lateral PSF at z = ', num2str(Media.MP(4, 3)), ' wavelengths'))
    % xlim([1, xnumpixv])
    % xlabel('x pixel')
    % ylabel('Intensity')
    % % 
    % figure
    % plot(abs_IQv_sum(:, 64, 1), 'LineWidth', 1.5)
    % title(strcat('Axial PSF at x = ', num2str(64 - numElements/2), ' wavelengths'))
    % % title('Axial PSF')
    % % xlim([d_pad_size(1), znumpix])
    % xlim([1, znumpixv])
    % xlabel('z pixel')
    % ylabel('Intensity')
    % 
    % figure
    % plot(abs_IQv_sum(:, 104, 1), 'LineWidth', 1.5)
    % title(strcat('Axial PSF at x = ', num2str(104 - numElements/2), ' wavelengths'))
    % % title('Axial PSF')
    % % xlim([d_pad_size(1), znumpix])
    % xlim([1, znumpixv])
    % xlabel('z pixel')
    % ylabel('Intensity')
    
    %% Plot PSFs on same graph
    
    % figure
    % yyaxis left
    % plot(apt_sum_intensity(697, :, 1), 'b', 'LineWidth', 1)
    % ylabel('Intensity')
    % yyaxis right
    % plot(abs_IQv_sum(697, :, 1), 'r--', 'LineWidth', 1)
    % title(strcat('Lateral PSF at z = ', num2str(Media.MP(4, 3)), ' wavelengths'))
    % xlim([1, xnumpix])
    % xlabel('x pixel')
    % ylabel('Intensity')
    % legend('Self recon', 'Verasonics recon')
    % % 
    % figure
    % yyaxis left
    % plot(apt_sum_intensity(:, 64, 1), 'b', 'LineWidth', 1)
    % ylabel('Intensity')
    % yyaxis right
    % plot(abs_IQv_sum(:, 64, 1), 'r--', 'LineWidth', 1)
    % title(strcat('Axial PSF at x = ', num2str(64 - numElements/2), ' wavelengths'))
    % % title('Axial PSF')
    % % xlim([d_pad_size(1), znumpix])
    % xlim([1, znumpix])
    % xlabel('z pixel')
    % ylabel('Intensity')
    % legend('Self recon', 'Verasonics recon')
    % 
    % 
    % figure
    % yyaxis left
    % plot(apt_sum_intensity(:, 104, 1), 'b', 'LineWidth', 1)
    % ylabel('Intensity')
    % yyaxis right
    % plot(abs_IQv_sum(:, 104, 1), 'r--', 'LineWidth', 1)
    % title(strcat('Axial PSF at x = ', num2str(104 - numElements/2), ' wavelengths'))
    % % title('Axial PSF')
    % % xlim([d_pad_size(1), znumpix])
    % xlim([1, znumpix])
    % xlabel('z pixel')
    % ylabel('Intensity')
    % legend('Self recon', 'Verasonics recon')
    
    
    %%
    % aa = angles(1);
    % % aa = 0;
    % k = 1 / wl;
    % phase_analytical_angle1 = zeros(znumpix, xnumpix);
    % for zp = 1:znumpix
    %     z_level_pix = zp / znumpix * zsize;      % z value of zi in m
    % 
    %     for xp = 1:xnumpix
    %         x_level_pix = (xp - xnumpix/2) / xnumpix * xsize;  % x value of xi in m
    %         phase_analytical_angle1(zp, xp) = mod(2 * k * (zp * cos(aa) + (sign(aa) * x_end - xp) * sin(aa) + d_correction_1way), 2*pi);
    %     end
    % end
    % figure;
    % imagesc(phase_analytical_angle1)
    %%
    %     else
    %     end
    
    
    % end

end
