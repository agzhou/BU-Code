function snrGG_final = calActivatedImage(mtrialsumGG, trial)

nthres1 = 5;
nthres2 = 10;
mskThreshold = 10;

snrGGnum1 = mean(mtrialsumGG(:,:,trial.nRest+2:trial.nRest+trial.nStim+1),3);
snrGGnum2 = mean(cat(3,mtrialsumGG(:,:,1:5),mtrialsumGG(:,:,end-15: end)),3);
snrGGnum = snrGGnum1 - snrGGnum2;
snrGGden = sqrt(5/25*var(mtrialsumGG(:,:,trial.nRest+1:trial.nRest+trial.nStim+1),1,3)+20/25*var(cat(3,mtrialsumGG(:,:,1:5),mtrialsumGG(:,:,end-15: end)), 1, 3));
% snrGGnum = mean(mtrialsumGG,3);
% snrGGden = std(mtrialsumGG, 1, 3);    
snrGG = snrGGnum1>snrGGnum2+snrGGden;

snrGG_final=bwareafilt(logical(snrGG),[mskThreshold,inf]); % 9 pixels connectivity filter 

%snrGG_final = snrGG_conv; % snrGGdB

end



function af=medfilter(a,n)
[nz,nx]=size(a);
af=a;
for iz=1+n:nz-n
    for ix=1+n:nx-n
        tmp=a([-n:n]+iz,[-n:n]+ix);
        st=std(tmp(:));
        m=median(tmp(:));   
        if abs(af(iz,ix)-m)>1.5*st
        af(iz,ix)=m;
        end
    end
end
end

