%% Description

% Compress unstacked row-column array (RC15gV) RF data to the receive channels the RCA beamformers use.
% Each transmit only uses half of the receive channels: the CR transmits (1:na, TX columns) receive on the rows and the RC
% transmits (na+1:2na, TX rows) receive on the columns (see grp1/grp2 in RCA_DAS_rcabeam, RCA_DAS, RCA_DAS_ffdas,
% RCA_DAS_fixed). The other half is dropped, and the kept channels are put in element order (P.Trans.Connector applied).

% Input:
%   r   unstacked RF, (# samples per frame, # receive channels, # frames), raw channel order, as from unstackFrames
%   P   acquisition parameter struct (Trans, Receive, na, numElements)
% Output:
%   rc  (# samples per acquisition, 2*na, numElements, # frames), same class as r (int16 from FastSave)
%         dim 2:  1:na      CR (TX columns), dim 3 = row elements    (remapped channels numElements+1 : 2*numElements)
%                 na+1:2na  RC (TX rows),    dim 3 = column elements (remapped channels 1 : numElements)
%   layout  struct describing rc, to save with it

% Notes:
%   - Each acquisition is cut out with its own Receive(a).startSample : endSample, so gaps between acquisitions in the
%     frame (if any) are dropped too
%   - Lossy for anything that reads the other channels (e.g., the Verasonics VSX recon in RCA_VerasonicsRecon_loop*.m, which
%     needs the full buffer layout), lossless for the RCA beamformers above
%   - RCA_DAS_rcabeam reads this layout directly (it detects it from the size)

function [rc, layout] = compressRCA_RF(r, P)

    N = P.numElements;
    na = P.na;
    Receive = P.Receive;
    if numel(Receive) < 2*na
        error('P.Receive has %d entries, need at least 2*na = %d', numel(Receive), 2*na)
    end
    ns = Receive(1).endSample - Receive(1).startSample + 1;
    for a = 2 : 2*na
        if Receive(a).endSample - Receive(a).startSample + 1 ~= ns
            error('Acquisition %d has %d samples, acquisition 1 has %d', a, Receive(a).endSample - Receive(a).startSample + 1, ns)
        end
    end
    if Receive(2*na).endSample > size(r, 1)
        error('The RF has %d samples per frame, but Receive(%d) ends at sample %d. Unstack the frames first (unstackFrames)', ...
            size(r, 1), 2*na, Receive(2*na).endSample)
    end
    conn = P.Trans.Connector(:)';
    if numel(conn) ~= 2*N || max(conn) > size(r, 2)
        error('P.Trans.Connector does not fit %d elements and %d receive channels', 2*N, size(r, 2))
    end

    chCR = conn(N+1 : 2*N); % raw channels of the row elements (RX of the CR transmits)
    chRC = conn(1 : N);     % raw channels of the column elements (RX of the RC transmits)

    nf = size(r, 3);
    rc = zeros(ns, 2*na, N, nf, 'like', r);
    for a = 1 : 2*na
        if a <= na, ch = chCR; else, ch = chRC; end
        rc(:, a, :, :) = reshape(r(Receive(a).startSample : Receive(a).endSample, ch, :), ns, 1, N, nf);
    end

    layout.format = 'RCA_compressed_v1';
    layout.dims = {'sample', 'acquisition', 'element', 'frame'};
    layout.description = ['(samples per acquisition, 2*na, numElements, frames). Acquisitions 1:na (CR) keep the row ' ...
        'elements, na+1:2*na (RC) the column elements, in element order (Trans.Connector already applied)'];
    layout.ns = ns;
    layout.na = na;
    layout.numElements = N;
    layout.numFrames = nf;

end
