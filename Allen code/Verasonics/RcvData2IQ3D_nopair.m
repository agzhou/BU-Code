%% Description

% For use with the 0 deg R/C combo plane wave tests
    
% Use delay-and-sum on RF data collected from a row-column array probe

% Input: RF data and various parameters
% Output: volumetric IQ data, phase, intensity, etc.

% Notes: 
%           - this uses padding for a nonzero start depth
%           - can use a gain for increasing depths according to some
%             known/estimated attenuation coefficient

function [IQ] = RcvData2IQ3D_nopair(RcvData, P)
    
    
    %%
    useGain = 1;
    % Currently it's set to use the element sensitivity as a weight
    
    %% potential to do/issues
    
    % I feel like I'm still missing something that improves quality
    
    % add band pass filter?
    
    %% Use parallel processing for speed
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
    numSubFrames = [];
    numSupFrames = [];

    assignStructVars(P)

    nf = numSubFrames;
    tw_peak = TW.peak;
    
%     assignFromParameterStructure;
%     wl
%     parfeval(@assignStructVars, 0, P); % numel(fieldnames(P))
    %% Get raw data and remap according to element-channel pairings
    tstart = clock;
    
    if iscell(RcvData)
        r = RcvData{1, 1}; % use this for the simulation where data is saved together at the end
    else
        r = RcvData;         % use this for continuous acquisition (ctsacq) when data is saved after each set of frames
    end
    
    %%%%%%%%%%%%%%%%
    % nf = 50;
    % r = r(:, :, 1:nf);
    %%%%%%%%%%%%%%%%
    
    A = Media.attenuation; % dB/cm/MHz
    A_per_m = A * Trans.frequency * 100; % dB/m
    
%     elemToConMap = Trans.Connector; % map channels to the transducer elements they're connected to
    c = Resource.Parameters.speedOfSound;
    
    r = r(:, Trans.Connector, :); % Reorganize data according to array element-to-channel map
%     clear elemToConMap % save memory
    s = size(r);        % (# z samples * # angles, # receive channels, # frames)
    
    %     t_correction = 0; % Transducer lens correction
    t_correction = 2* (Trans.lensCorrection * wl) / c; % Trans.lensCorrection is the one way delay in wavelengths through lens
    d_correction_1way = (Trans.lensCorrection * wl);     % use further samples to account for delay.......
    
    element_sens_cutoff = 0.6; %%%%%%%%%%%%%%%%%%%%%%%%%
    samplesPerWL = Receive(1).samplesPerWave;
    timePerSample = 1/(Trans.frequency * 1e6)/samplesPerWL;
    distPerSample = wl/samplesPerWL;
    
    %     if length(s) == 3 % 1D array/2D data
    
    %% reorganize RcvData
    
    % Don't need the full 3D space yet, since value of rcv voltage at each
    % element is same along the length
    nzs_nopad = s(1) / na / nf;
    nxs_nopad = numElements;
    nys_nopad = nxs_nopad;
    d = zeros(nzs_nopad, numElements, na, nf); % reorganize RcvData into (# z samples, # x or y samples, # acquisitions (# angles), # frames)

    for f = 1:nf
    
        switch f
%             case 1
%                 % column rcv
%                 d(:, :, 1, f)  = r(Receive(f).startSample:Receive(f).endSample, numElements + 1 : numElements * 2);
%             case 2
%                 % row rcv
%                 d(:, :, 1, f) = r(Receive(f).startSample:Receive(f).endSample, 1:numElements);
%             case 3
%                 % row rcv
%                 d(:, :, 1, f) = r(Receive(f).startSample:Receive(f).endSample, 1:numElements);
%             case 4
%                 % column rcv
%                 d(:, :, 1, f)  = r(Receive(f).startSample:Receive(f).endSample, numElements + 1 : numElements * 2);

            case 1
                % column rcv
                d(:, :, 1, f)  = r(Receive(f).startSample:Receive(f).endSample, 1:numElements);
            case 2
                % row rcv
                d(:, :, 1, f) = r(Receive(f).startSample:Receive(f).endSample, numElements + 1 : numElements * 2);
            case 3
                % row rcv
                d(:, :, 1, f) = r(Receive(f).startSample:Receive(f).endSample, numElements + 1 : numElements * 2);
            case 4
                % column rcv
                d(:, :, 1, f)  = r(Receive(f).startSample:Receive(f).endSample, 1:numElements);
        end
    
    end
%     d = squeeze(d);

    % (# x samples, # y samples, # z samples, # acquisitions (# angles), # frames)
    d = permute(d, [2, 1, 3, 4]); 
    
    % d_old = d; %%%%% for testing purposes
    d_pad = zeros(nxs_nopad, round(startDepth * samplesPerWL), na, nf);
    d_pad_size = size(d_pad);
    d = [d_pad, d, d_pad];
    
    sd = size(d);
    nxs = sd(1);
    nys = sd(2);
    nzs = sd(3);
    
    clear r % save memory

%     figure; imagesc(squeeze(d(:, :, 1, 1))')
    
    %         figure
    %         imagesc(squeeze(d(:, :, 1, 1)))
    %         figure
    %         imagesc(squeeze(d(:, 80, :, 1, 1)))
    %         figure
    %         imagesc(squeeze(d(:, :, 1, 2, 1)))
    
    %% Plot slices (for checking)
    % xp = sd(1);
    % yp = sd(2);
    % zp = sd(3);
    % % [X, Y, Z] = meshgrid(1:xp, 1:yp, 1:zp);
    % % [X, Y, Z] = meshgrid(1:zp, linspace(-yp/2, yp/2, yp), linspace(-xp/2, xp/2, xp));
    % [X, Y, Z] = meshgrid(linspace(-xp/2, xp/2, xp), linspace(-yp/2, yp/2, yp), 1:zp);
    % 
    % % X = permute(X, [2, 1, 3]);
    % % Y = permute(Y, [2, 1, 3]);
    % % Z = permute(Z, [2, 1, 3]);
    % 
    % % X = permute(X, [1, 3, 2]);
    % % Y = permute(Y, [1, 3, 2]);
    % % Z = permute(Z, [1, 3, 2]);
    % % [Z, X, Y] = ndgrid(1:zp, linspace(-xp/2, xp/2, xp), linspace(-yp/2, yp/2, yp));
    % 
    % % ptZ = unique(Media.MP(:, 3)); % all z points of scatterers from simulation
    % 
    % ptZ = 50;
    % % xslice = xp/2 + 1;   % middle
    % % yslice = yp/2 + 1;   % middle
    % xslice = 0;   % middle
    % yslice = 0;   % middle
    % % zslice = ptZ*2;
    % zslice = ceil((ptZ ./ endDepth) .* zp);
    % 
    % ddd = squeeze(d(:, :, :, 1, 1));
    % % ddd = permute(ddd, [2, 3, 1]);
    % % size(ddd)
    % figure;
    % % slice(X, Y, Z, ddd, xslice, yslice, zslice)
    % % slice(Z, X, Y, ddd, xslice, yslice, zslice)
    % slice(X, Y, Z, ddd, xslice, yslice, zslice)
    % 
    % shading flat
    % title(['Original frame ' num2str(1)])
    % xlabel('x pixels')
    % ylabel('y pixels')
    % zlabel('z pixels')
    
    
    
    
    %         imgv = ImgData{1, 1};
    %         imagesc(squeeze(imgv(:, 1, :, 2)))
    %         imagesc(squeeze(imgv(:, :, 50, 1)))
    
    
    %%%%%%%%%%%%%%%%%%%%
    
            % reorganize RcvData into (# z samples, # angles, # x samples, # frames)
            % angles is second bc that's how Verasonics stacks it
    %         d = reshape(r, s(1)/na, na, s(2), s(3));
    % %         dtest = reshape(d, s(1)/na, s(2), na, s(3));
    
    %         figure
    %         imagesc(squeeze(d(:, 1, :, 1)))
    %         figure
    %         imagesc(squeeze(d(:, 2, :, 1)))
    
    
    %         distTotal = repmat(distTotal_single', 1, Trans.numelements);
    
    %% Define reconstruction region and get time delays for each pixel
    pixsizez = wl/1;
    pixsizex = wl;
    pixsizey = wl;
    
    xsize = numElements * Trans.spacingMm / 1e3; % region's x in m
    ysize = numElements * Trans.spacingMm / 1e3; % region's y in m
    zsize = (endDepth - startDepth) * wl;    % region's z in m
    
    % znumpix = ceil(zsize / pixsizez);         % # rows
    % xnumpix = ceil(xsize / pixsizex);         % # cols
    % ynumpix = ceil(ysize / pixsizey);         % # cols
    
    xnumpix = nxs_nopad*1;
    ynumpix = nys_nopad*1;
    znumpix = nzs_nopad;
    
    arrayLength = numElements * Trans.spacingMm / 1e3;
    % elementx = - (arrayLength - (Trans.spacingMm / 1e3)) / 2 : Trans.spacingMm / 1e3 : (arrayLength - (Trans.spacingMm / 1e3)) / 2; % define x positions of array elements in m, with x = 0 at center
    % elementx = linspace((1 - numElements/2 - 1/2) / numElements * xsize, (numElements - numElements/2 - 1/2) / numElements * xsize, numElements);
    % elementx = linspace((1 - numElements/2 - 1/2) / numElements * xsize, (numElements - numElements/2 - 1/2) / numElements * xsize, numElements);
    % elementx = linspace((1 - numElements/2) / numElements * xsize, (numElements - numElements/2) / numElements * xsize, numElements);
    
    % elementx = - (arrayLength - (Trans.spacingMm / 1e3)) / 2 : Trans.spacingMm / 1e3 : (arrayLength - (Trans.spacingMm / 1e3)) / 2; % define x positions of array elements in m, with x = 0 at center
    % elementy = elementx;
    elementx = (wl .* Trans.ElementPos(1:numElements, 1))';
    elementy = (wl .* Trans.ElementPos(numElements + 1 : numElements*2, 2))';
    x_end = max(elementx);
    y_end = max(elementy);
    
    % distTotal_single = distPerSample .* (1 : s(1)/na / 2); % z distance in m, at every z sample point
    
    d_sum = zeros(xnumpix, ynumpix, znumpix, na, nf, 'single'); %%%%%%%%%%%%%%
    %                 d_sum = zeros(znumpix, na, xnumpix, nf);
    
    d_delayed_temp = zeros(numElements, nf); % get the value at a single sample point, supposed to be the max of the envelope
                                       % (# x or y samples, # frames)
    
    % m = double(max(distTotal_single));
    
    %%%% make this more efficient?
    
    frame_ind_temp = zeros(numElements * nf * na, 1);
    for f = 1:nf
        frame_ind_temp((f-1)*numElements * na + 1 : f * na * numElements) = f;
    end
    
    angles_ind_temp = zeros(numElements * na, 1);
    % first half of angles
    for a = 1:na
        angles_ind_temp((a-1)*numElements + 1 : a * numElements) = a;
    end
    
    angles_ind_temp = repmat(angles_ind_temp, nf, 1);
    

    %% Get travel delays
    
    % Go through each pixel (xp, yp, zp)
    parfor zp = 1:znumpix % do z first so parfor is more efficient, theoretically
%     for zp = 600
%         disp('test')
%         disp(wl)
%         disp(zp)
    % for zp = round(70/endDepth*znumpix)
%         assignStructVars(P)
    %             z_level_pix = zp / znumpix * zsize;      % z value of zp in m. z = 0 is at the level of the transducer array
        z_level_pix = startDepth * wl + (zp - 1/2) / znumpix * zsize;      % z value of zp in m. z = 0 is at the level of the transducer array

        for yp = 1:ynumpix
%         for yp = 20
    %     for yp = Media.MP(1, 2) + numElements/2
    %         y_level_pix = (yp + ynumpix/2 - 1/2) / ynumpix * ysize;  % y value of yp in m. y = 0 is center of transducer array
    
            y_level_pix = (yp - (ynumpix + 1)/2) / ynumpix * ysize;  % y value of yp in m. y = 0 is center of transducer array
    
            for xp = 1:xnumpix
%             for xp = 40
    %         for xp = Media.MP(1, 1) + numElements/2 
            % xp counts positively from -numElements/2

                x_level_pix = (xp - (xnumpix + 1)/2) / xnumpix * xsize;  % x value of xp in m. x = 0 is center of transducer array
    
                % RX travel dist from pixel to each transducer element e
                % dimensions (# (y or x) elements, # acquisitions)
    
                % when TX angle rotates about y axis, it's the opposite (x) elements that
                % receive, so calculate distance with the y coord for
                % traveldistRX_y_temp
    %             traveldistRX_temp = sqrt(z_level_pix ^ 2 + repmat((elementx - x_level_pix).^2, numElements, 1) + repmat((elementy - y_level_pix).^2, numElements, 1)');
                traveldistRX_y_temp = (sqrt(z_level_pix ^ 2 + (elementy - y_level_pix).^2) + 1 * d_correction_1way)';
                traveldistRX_x_temp = (sqrt(z_level_pix ^ 2 + (elementx - x_level_pix).^2) + 1 * d_correction_1way)';
                
                % (# elements, # frames)
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 traveldistRX_temp = [traveldistRX_y_temp, traveldistRX_x_temp, traveldistRX_x_temp, traveldistRX_y_temp];
                traveldistRX_temp = [traveldistRX_x_temp, traveldistRX_y_temp, traveldistRX_y_temp, traveldistRX_x_temp];
    
                % travel dist from TX to pixel depends on z and angle
                % traveldistTX_temp has dimensions (numElements, # acquisitions) 
                % where the values in each column are the same
    %             traveldistTX_temp_y = z_level_pix .* cos(angles) + (sign(angles) .* y_end - y_level_pix) .* sin(angles) + 1 * d_correction_1way;    % First half of acqs are for angles changing in y
    %             traveldistTX_temp_x = z_level_pix .* cos(angles) + (sign(angles) .* x_end - x_level_pix) .* sin(angles) + 1 * d_correction_1way;    % Second half of acqs are for angles changing in x
       
                % when angle is rotating about y axis, it's the x elements
                % doing that, so use x for traveldistTX_temp_y
                traveldistTX_temp_y = z_level_pix .* cos(angles) + (x_level_pix + sign(angles) .* x_end) .* sin(angles) + 1 * d_correction_1way;    % First half of acqs are for angles changing in y
                traveldistTX_temp_x = z_level_pix .* cos(angles) + (sign(angles) .* y_end - y_level_pix) .* sin(angles) + 1 * d_correction_1way;    % Second half of acqs are for angles changing in x
                traveldistTX_temp_x = fliplr(traveldistTX_temp_x);

                traveldistTX_temp = repmat([traveldistTX_temp_y, traveldistTX_temp_x, traveldistTX_temp_x, traveldistTX_temp_y], numElements, 1);
%                 traveldistTX_temp = repmat([traveldistTX_temp_x, traveldistTX_temp_y, traveldistTX_temp_y, traveldistTX_temp_x], numElements, 1);
                
                % (# x elements, # y elements, # acquisitions)
    %             temp_td = traveldistTX_temp_matrix + repmat(traveldistRX_temp, 1, 1, na*2);
    %             temp_td = traveldistTX_temp + traveldistRX_temp + wl* (TW.Numpulses/2);
                temp_td = traveldistTX_temp + traveldistRX_temp + wl* (tw_peak); % TW.peak is the # of wavelengths to peak of transmitted waveform

                %%%%%%%%%%%%%%%%%%%%%%%%%%%% wtf %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 temp_td = temp_td(:, 1);
    
                % do I want to unaccount for the angle if the pixel is out of the beam
    
                % ind has dimensions (# (y or x) elements, # acquisitions)
    %             ind = round(temp_td ./ m .* length(distTotal_single)); % get z sample indices for envelope values corresponding to each element e, at a single pixel
                ind = round(temp_td ./ distPerSample - startDepth*samplesPerWL);
    
%                 ia = repmat(reshape(ind, numElements*na, 1), nf, 1); %%%% rearrange indices
%                 ia = repmat(ind, nf, 1); %%%% rearrange indices
                ia = ind(:);
                ib = repmat((1:numElements)', nf*na, 1);
    
                ind_linear = sub2ind(sd, ib, ia, angles_ind_temp, frame_ind_temp);
    
                d_delayed_temp = d(ind_linear); % vector that stacks the delayed d value from each frame            
    
                % unstack the delayed d values
%                 dd = zeros(nxs, na, nf); 
                dd = reshape(d_delayed_temp, numElements, numSubFrames);
%                 a = 1;
%                 for f = 1:nf %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     % First half of angles (rotate about y)
%                     
%                         dd(:, a, f) = d_delayed_temp((f-1)*(na) * numElements + (a-1)*numElements + 1 : (f-1)*(na)*numElements + (a)*numElements);
%         
%                     
% %                         dd(:, na + a, f) = d_delayed_temp((f-1)*(na * 2) * numElements + (na + a-1)*numElements + 1 : (f-1)*(na * 2)*numElements + (na + a)*numElements);
% 
%                 end
        
%                 dd = squeeze(dd);
    
                % Account for element angle sensitivity 
                % angle is from a pixel to each element
                angle_p_y = atan((y_level_pix - elementy) ./ z_level_pix); % for acq where angle is rotating about y, which means RX has the x elements active
                angle_p_x = atan((x_level_pix - elementx) ./ z_level_pix);
                % interpolate provided element vs. angle sensitivity curve
    %             element_sens_p_y = interp1(linspace(-pi/2, pi/2, length(Trans.ElementSens)), Trans.ElementSens, angle_p_y);
    %             element_sens_p_x = interp1(linspace(-pi/2, pi/2, length(Trans.ElementSens)), Trans.ElementSens, angle_p_x);
    
                element_sens_p_y = interp1(linspace(-pi/2, pi/2, length(Trans.ElementSens)), Trans.ElementSens, angle_p_y)';
                element_sens_p_x = interp1(linspace(-pi/2, pi/2, length(Trans.ElementSens)), Trans.ElementSens, angle_p_x)';
                
                % commented this on 10/23/24
    %             element_sens_mask_y = element_sens_p_y > element_sens_cutoff;
    %             element_sens_mask_x = element_sens_p_x > element_sens_cutoff;
    %             dd([repmat(~element_sens_mask_y', 1, na, nf), repmat(~element_sens_mask_x', 1, na, nf)]) = 0;
    
                % new 10/23/24 11/7/24
%                 dd = dd .* [repmat(element_sens_p_y, 1, na, nf), repmat(element_sens_p_x, 1, na, nf)]; % element sensitivity weighting



%                 dd = dd .* [element_sens_p_y, element_sens_p_x, element_sens_p_x, element_sens_p_y]; % element sensitivity weighting
                
                if useGain
                    dd = dd .* 10.^((temp_td .* -A_per_m/1)  ./ 20);
                end
        
                d_sum(xp, yp, zp, :, :) = squeeze(sum(dd, 1));
            end
        end
        
    end
    
    % Do I need to transpose for y vs x???
    
    
    %%
    
    tend = clock;
    disp(etime(tend, tstart))
    
    % Calculate IQ
    d_sum = permute(d_sum, [3, 1, 2, 4, 5]); 
    
    % IQ = permute(hilbert(d_sum_p), [2, 3, 1, 4, 5]); % operates on the columns (z) of the d_sum_p matrix
    % change 9/25/24
    IQ = permute(hilbert(d_sum), [3, 2, 1, 4, 5]); % operates on the columns (z) of the d_sum_p matrix

    
end
