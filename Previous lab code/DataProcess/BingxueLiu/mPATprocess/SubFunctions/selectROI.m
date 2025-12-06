function bmsk = selectROI(X)
figure; imshow(X,[],'InitialMagnification','fit');
roi = drawfreehand('Smoothing', 10);
bmsk = createMask(roi);
end

