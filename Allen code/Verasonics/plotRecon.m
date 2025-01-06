%% changes

% 11/7/24: added linear array functionality




%%
function plotRecon(IQ, P, fn)

%     phase = atan2(imag(IQ), real(IQ));     % Need atan2 for full 0 - 2pi range!!
   
    samplesPerWL =  P.Receive(1).samplesPerWave;
    
    if numel(size(IQ)) == 4 % linear array
        %% Coherently sum across angles
        
        IQ_coherent_sum = squeeze(sum(IQ, 3));
        I_coherent_sum = abs(IQ_coherent_sum); % intensity

        %% Plot plane
        figure
        imagesc(squeeze(I_coherent_sum(:, :, fn)));

    else                    % row column array

        %% Coherently sum across angles
        
        IQ_coherent_sum = squeeze(sum(IQ, 4));
        I_coherent_sum = abs(IQ_coherent_sum); % intensity
       
        %% Plot planes
    %     fn = 3;
        figure
        imagesc(squeeze(I_coherent_sum(:, 40, :, fn))')
        title('yz plane at x = 0 mm')
        ylabel('z pixels')
        xlabel('y pixels')
        
        figure
        imagesc(squeeze(I_coherent_sum(40, :, :, fn))')
        title('xz plane at y = 0 mm')
        ylabel('z pixels')
        xlabel('x pixels')
        
        % figure
        % imagesc(squeeze(I_coherent_sum(:, 68, :, fn))')
        % title('yz plane at x = 3 mm')
        % ylabel('z pixels')
        % xlabel('y pixels')
        % figure
        % imagesc(squeeze(I_coherent_sum(40, :, :, fn))')
        
        % xy
        % figure
        % imagesc(squeeze(I_coherent_sum(:, :, 697, fn))')
        % title('xy plane at z = 10 mm')
        % ylabel('y pixels')
        % xlabel('x pixels')
        % axis square
        % 
        % figure
        % imagesc(squeeze(I_coherent_sum(:, :, 1392, fn))')
        % title('xy plane at z = 20 mm')
        % ylabel('y pixels')
        % xlabel('x pixels')
        % axis square
        % 
        % figure
        % imagesc(squeeze(I_coherent_sum(:, :, 2086, fn))')
        % title('xy plane at z = 30 mm')
        % ylabel('y pixels')
        % xlabel('x pixels')
        % axis square
        
        %% Plot volume with slices after reconstruction
    %     xpn = xnumpix;
    %     ypn = ynumpix;
    %     zpn = znumpix;
    %     % [X, Y, Z] = meshgrid(1:xp, 1:yp, 1:zp);
    %     [X, Y, Z] = meshgrid(linspace(-ypn/2, ypn/2, ypn), linspace(-xpn/2, xpn/2, xpn), 1:zpn);
    %     % [X, Y, Z] = meshgrid(linspace(-xpn/2, xpn/2, xpn), linspace(-ypn/2, ypn/2, ypn), 1:zpn);
    %     
    %     % [X, Y, Z] = meshgrid(linspace(xpn/2, -xpn/2, xpn), linspace(-ypn/2, ypn/2, ypn), 1:zpn);
    %     
    %     ptZ = unique(Media.MP(:, 3)); % all z points of scatterers from simulation
    %     
    %     
    %     % xslice = 0;   % middle
    %     % yslice = 0;   % middle
    %     
    %     xslice = unique(Media.MP(:, 1));
    %     yslice = unique(Media.MP(:, 2));
    %     zslice = ceil((ptZ ./ endDepth) .* zpn); %%%% fix this, for some reason it's not exactly on the points
    %     
    %     % plot frame 1 coherent sum
    %     figure;
    %     fn = 2;
    %     % testa = permute(I_coherent_sum, [2, 1, 3, 4, 5]);
    %     slice(X, Y, Z, squeeze(I_coherent_sum(:, :, :, fn)), xslice, yslice, zslice)
    %     % slice(X, Y, Z, squeeze(testa(:, :, :, fn)), xslice, yslice, zslice)
    %     
    %     shading flat
    %     title(['Self recon frame ' num2str(fn)])
    %     xlabel('x pixels')
    %     ylabel('y pixels')
    %     zlabel('z pixels')
    %     alpha('color');
    %     colormap(jet);
        
        
        %% Plot Verasonics recon from sim in a slice
        % IQv = squeeze(IData{1, 1} + 1i.*QData{1, 1});
        % abs_IQv_sum = squeeze(abs(sum(IQv, 4)));
        % 
        % [xpv, ypv, zpv, ~] = size(abs_IQv_sum);
        % [Xv, Yv, Zv] = meshgrid(linspace(-ypv/2, ypv/2, ypv), linspace(-xpv/2, xpv/2, xpv), 1:zpv);
        % 
        % 
        % %
        % zslicev = ceil((ptZ ./ endDepth) .* zpv); %%%% fix this, for some reason it's not exactly on the points
        % figure
        % slice(Xv, Yv, Zv, squeeze(abs_IQv_sum(:, :, :, fn)), xslice, yslice, zslicev)
        % 
        % shading flat
        % title(['Verasonics recon frame ' num2str(fn)])
        % xlabel('x pixels')
        % ylabel('y pixels')
        % zlabel('z pixels')      
        % alpha('color');
        % colormap(jet);
        
        %% verasonics slices
        
        % fn = 1;
        % figure
        % imagesc(squeeze(abs_IQv_sum(:, 40, :, fn))')
        % title('Verasonics yz plane at x = 0 mm')
        % ylabel('z pixels')
        % xlabel('y pixels')
        % 
        % figure
        % imagesc(squeeze(abs_IQv_sum(40, :, :, fn))')
        % title('Verasonics xz plane at y = 0 mm')
        % ylabel('z pixels')
        % xlabel('x pixels')
        %% Volume Viewer
        % fn = 1;
        volumeViewer(squeeze(I_coherent_sum(:, :, :, fn)), scaleFactors = [1, 1, 1/samplesPerWL])
        
        % IQv = squeeze(IData{1, 1} + 1i.*QData{1, 1});
        % abs_IQv_sum = squeeze(abs(sum(IQv, 4)));
        
        % IQv1f = squeeze(IQv(:, :, :, :, 1));
        % abs_IQv1f_sum = squeeze(abs(sum(IQv1f, 4)));
        
        % figure
        % imagesc(squeeze(abs_IQv_sum(:, 40, :, fn))')
        % imagesc(squeeze(abs_IQv_sum(Media.MP(2, 1) + numElements/2 , :, :, fn))')
        % volumeViewer(squeeze(abs_IQv1f_sum(:, :, :, fn)), scaleFactors = [1, 1, 1/samplesPerWL])
        
        %%
        % fn = 1;
        % figure
        % imagesc(squeeze(abs_IQv_sum(:, 40, :, fn))')
        % figure
        % imagesc(squeeze(abs_IQv_sum(40, :, :, fn))')
    end
end