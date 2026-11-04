function P = cwtPrepare(fb, cls)
% CWTPREPARE  Pull the frequency-domain filters and padding out of a cwtfilterbank (once).
%   P.psiT  : [nPad x nFreq] filters, in class cls
%   P.pad   : samples of boundary extension on each side
%   P.reflect : true for 'reflection' boundary, false for 'zeros'
%   P.faxis : [nFreq x 1] center frequencies [Hz]
% Uses struct(fb) because PsiDFT/SignalPad are hidden properties (may change between releases);
% the result is validated against cwt() on a test signal and errors out if they disagree.
nT = fb.SignalLength;
t  = cast(sin(2*pi*(1:nT).'/37) + (1:nT).'/nT, cls);
[wref, faxis] = cwt(t, 'FilterBank', fb);          % also primes fb.PsiDFT in class cls
w = warning('off', 'all'); st = struct(fb); warning(w);
P.psiT    = cast(st.PsiDFT, cls).';
P.pad     = st.SignalPad;
P.reflect = startsWith(fb.Boundary, 'r');
P.faxis   = faxis;
% self-check
c = cwtBatchedCore(t, P);                           % [nT x 1 x nF]
c = permute(c, [3 1 2]);
err = max(abs(c - wref), [], 'all');
if err > 1e-4 * max(abs(wref), [], 'all')
    error('cwtPrepare:mismatch', 'Manual CWT disagrees with cwt() (err %g); fall back to cwt().', err);
end
end
