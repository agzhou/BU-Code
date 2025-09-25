% Description:
%       Overlay 3D activation maps onto a template (like CBV index)

%%
compareVolumes(CBViallSFAvg, am_rPDI)
% compareVolumes(CBViallSFAvg, am_rCBV)
% compareVolumes(CBViallSFAvg, am_rCBFspeed)








%%
generateTiffStack_multi({CBVi_allSF_avg .^ 0.7}, [8.8, 8.8, 8], 'gray', 1)


vcmap = colormap_ULM;

generateTiffStack_multi({am_rPDI}, [8.8, 8.8, 8], 'jet', 1)
% generateTiffStack_multi({am_rCBV}, [8.8, 8.8, 8], 'jet', 1)
generateTiffStack_multi({am_rCBV}, [8.8, 8.8, 8], vcmap, 1)
generateTiffStack_multi({am_rCBFspeed}, [8.8, 8.8, 8], 'jet', 1)
%%
function compareVolumes(vol1, vol2) % Can change this so it has a cell array input and goes through more than 2 volumes

    viewerThresholded = viewer3d(BackgroundColor = "black", BackgroundGradient="off");
    volshow(vol1 .^ 1, Parent=viewerThresholded, RenderingStyle = "MaximumIntensityProjection", ...
        Colormap=[1 0 1],Alphamap = "linear");
    volshow(vol2 .^ 1, Parent=viewerThresholded, RenderingStyle = "MaximumIntensityProjection", ...
        Colormap=[0 1 0],Alphamap = "linear");

end