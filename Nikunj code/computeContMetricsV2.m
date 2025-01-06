function [C, CNR, GCNR] = computeContMetricsV2(images,signalMap,noiseMap,imgMap)

% Input: 
%   
%   images - a struct array containing two fields. The first field is the
%   name of the corresponding image. The second field is the acutal image
%   data.
%   signalMap - logical array of 1s and 0s indicating signal subregion
%   noiseMap - logical array of 1s and 0s indicating noise subregion
%   imgMap - logical array of 1s and 0s indicating image subregion

% Output:
%   C - Contrast
%   CNR - Contrast to Noise Ratio
%   GCNR - Generalized Contrast to Noise Ratio


    nImg = length(images);
    C = zeros(nImg,1); CNR = C; GCNR = C;

    for n = 1:nImg
        img = images(n).data;

        mu_i=mean(img(signalMap)); mu_o=mean(img(noiseMap));
        v_i=var(img(signalMap)); v_o=var(img(noiseMap));

        C(n)=10*log10(mu_i./mu_o);
        CNR(n)=10*log10(abs(mu_i-mu_o)/sqrt(v_i+v_o));

        % Select subset of image to generate histogram for. Compute w.r.t.
        % magnitudes
        imgReg = img(imgMap);
        x=linspace(min(imgReg(:)),max(imgReg(:)),25);


        % [~,edges] = histcounts(img(signalMap), 'Normalization','pdf');
        % pdf_i = edges(2:end) - (edges(2)-edges(1))/2;
        % 
        % [~,edges] = histcounts(img(noiseMap), 'Normalization','pdf');
        % pdf_o = edges(2:end) - (edges(2)-edges(1))/2;

        [pdf_i]=hist(img(signalMap),x);
        [pdf_o]=hist(img(noiseMap),x);

        OVL=sum(min([pdf_i./sum(pdf_i); pdf_o./sum(pdf_o)]));
        GCNR(n)= 1 - OVL;

    end

end