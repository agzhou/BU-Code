% No frequency separation. Same inputs/outputs as IQ2g1T_3D_nofreqsep (current signature: IQ is already cropped by the
% caller, and there is no unused noise output), but the SVD clutter filter and the high pass filter are done together by
% clutterHPF_fast. g1T_fft is unchanged.

function [PDI, g1] = IQ2g1T_3D_nofreqsep_fast(IQ, sv_threshold_lower, sv_threshold_upper, HPF, nTau)

    IQf_HPF = clutterHPF_fast(IQ, sv_threshold_lower, sv_threshold_upper, HPF);
    clearvars IQ

    PDI = sum(abs(IQf_HPF).^2, 4) ./ size(IQf_HPF, 4);
    g1 = g1T_fft(IQf_HPF, nTau);
end
