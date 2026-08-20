function [vesselness, angleDeg] = vesselnessAngle2D(I, sigmas, spacing, tau, brightondark)
% calculates the vesselness probability map (local tubularity) of a 2D
% input image, together with the local vessel orientation (tangent
% direction) at the scale that produced each pixel's vesselness response.
%
% [vesselness, angleDeg] = vesselnessAngle2D(I, sigmas, spacing, tau, brightondark)
%
% This is a fork of vesselness2D.m (T. Jerman, 2014; based on D. Kroon's
% Frangi filter, 2009) that additionally tracks, at each scale, the
% eigenvector of the SMALLER-magnitude Hessian eigenvalue (Lambda1) --
% i.e. the direction of minimal curvature, which runs along a tubular
% structure -- the same way Kroon's eig2image.m does. Whichever scale
% wins the max-response vesselness at a pixel also supplies that pixel's
% angle, so vesselness and angleDeg are always scale-consistent.
%
% inputs,
%   I : 2D image
%   sigmas : vector of scales on which the vesselness is computed
%   spacing : input image spacing resolution [dim1; dim2] - during hessian
%       matrix computation, the gaussian filter kernel size in each
%       dimension can be adjusted to account for different image spacing
%       for different dimensions
%   tau : (between 0.5 and 1) : parameter that controls response uniformity
%       - lower tau -> more intense output response
%   brightondark: (true/false) : are vessels (tubular structures) bright on
%       dark background or dark on bright (default for 2D is false)
%
% outputs,
%   vesselness: maximum vesselness response over scales sigmas, in [0, 1]
%   angleDeg: local vessel (tangent) orientation [degrees], wrapped to
%       [-90, 90) since a line's orientation is undirected (mod 180 deg).
%       Measured as atan2d(dim2 component, dim1 component) of the tangent
%       eigenvector, i.e. 0 deg = running along dim 1 (rows), +-90 deg =
%       running along dim 2 (columns). For vUS_2D.m's PDI image
%       ([z voxels, x voxels]), that means 0 deg = along the depth/beam
%       (z) axis and +-90 deg = purely lateral (x). Pixels with
%       vesselness ~ 0 are set to NaN.
%
% example:
%   [V, A] = vesselnessAngle2D(I, 1:5, [1;1], 1, false);
%
% Function written by T. Jerman, University of Ljubljana (October 2014)
% Based on code by D. Kroon, University of Twente (May 2009)
% Orientation output added (adapted from D. Kroon's eig2image.m)

verbose = 1;

if nargin<5
    brightondark = false; % default mode for 2D is dark vessels compared to the background
end

I = single(I);

bestIx = zeros(size(I), 'single');
bestIy = zeros(size(I), 'single');

for j = 1:length(sigmas)

    if verbose
        disp(['Current filter scale (sigma): ' num2str(sigmas(j)) ]);
    end

    [~, Lambda2, Ix, Iy] = imageEigenvalues(I,sigmas(j),spacing,brightondark);
    if brightondark == true
        Lambda2 = -Lambda2;
    end

    % proposed filter at current scale
    Lambda3 = Lambda2;

    Lambda_rho = Lambda3;
    Lambda_rho(Lambda3 > 0 & Lambda3 <= tau .* max(Lambda3(:))) = tau .* max(Lambda3(:));
    Lambda_rho(Lambda3 <= 0) = 0;
    response = Lambda2.*Lambda2.*(Lambda_rho-Lambda2).* 27 ./ (Lambda2 + Lambda_rho).^3;

    response(Lambda2 >= Lambda_rho./2 & Lambda_rho > 0) = 1;
    response(Lambda2 <= 0 | Lambda_rho <= 0) = 0;
    response(~isfinite(response)) = 0;

    %max response over multiple scales
    if(j==1)
        vesselness = response;
        bestIx = Ix;
        bestIy = Iy;
    else
        upd = response > vesselness;
        vesselness(upd) = response(upd);
        bestIx(upd) = Ix(upd);
        bestIy(upd) = Iy(upd);
    end

    clear response Lambda2 Lambda3 Ix Iy
end

vesselness = vesselness ./ max(vesselness(:)); % should not be really needed
vesselness(vesselness < 1e-2) = 0;

angleDeg = mod(atan2d(bestIy, bestIx) + 90, 180) - 90; % wrap to [-90, 90)
angleDeg(vesselness == 0) = NaN;

function [Lambda1, Lambda2, Ix, Iy] = imageEigenvalues(I,sigma,spacing,brightondark)
% calculates the two eigenvalues (and the Lambda1 eigenvector, i.e. the
% tangent direction) for each pixel in an image

% Calculate the 2D hessian
[Hxx, Hyy, Hxy] = Hessian2D(I,sigma,spacing);

% Correct for scaling
c=sigma.^2;
Hxx = c*Hxx;
Hxy = c*Hxy;
Hyy = c*Hyy;

% reduce computation by computing vesselness only where needed
% S.-F. Yang and C.-H. Cheng, "Fast computation of Hessian-based
% enhancement filters for medical images," Comput. Meth. Prog. Bio., vol.
% 116, no. 3, pp. 215-225, 2014.
B1 = - (Hxx+Hyy);
B2 = Hxx .* Hyy - Hxy.^2;

T = ones(size(B1));

if brightondark == true
    T(B1<0) = 0;
    T(B2==0 & B1 == 0) = 0;
else
    T(B1>0) = 0;
    T(B2==0 & B1 == 0) = 0;
end

clear B1 B2;

indeces = find(T==1);

Hxx = Hxx(indeces);
Hyy = Hyy(indeces);
Hxy = Hxy(indeces);

% Calculate eigen values and the Lambda1 eigenvector (tangent direction)
[Lambda1i,Lambda2i,Ixi,Iyi]=eigvalOfHessian2D(Hxx,Hxy,Hyy);

clear Hxx Hyy Hxy;

Lambda1 = zeros(size(T));
Lambda2 = zeros(size(T));
Ix = zeros(size(T));
Iy = zeros(size(T));

Lambda1(indeces) = Lambda1i;
Lambda2(indeces) = Lambda2i;
Ix(indeces) = Ixi;
Iy(indeces) = Iyi;

% some noise removal
Lambda1(~isfinite(Lambda1)) = 0;
Lambda2(~isfinite(Lambda2)) = 0;

Lambda1(abs(Lambda1) < 1e-4) = 0;
Lambda2(abs(Lambda2) < 1e-4) = 0;


function [Dxx, Dyy, Dxy] = Hessian2D(I,Sigma,spacing)
%  filters the image with an Gaussian kernel
%  followed by calculation of 2nd order gradients, which aprroximates the
%  2nd order derivatives of the image.
%
% [Dxx, Dyy, Dxy] = Hessian2D(I,Sigma,spacing)
%
% inputs,
%   I : The image, class preferable double or single
%   Sigma : The sigma of the gaussian kernel used. If sigma is zero
%           no gaussian filtering.
%   spacing : input image spacing
%
% outputs,
%   Dxx, Dyy, Dxy: The 2nd derivatives

if nargin < 3, Sigma = 1; end

if(Sigma>0)
    F=imgaussian(I,Sigma,spacing);
else
    F=I;
end

% Create first and second order diferentiations
Dy=gradient2(F,'y');
Dyy=(gradient2(Dy,'y'));
clear Dy;

Dx=gradient2(F,'x');
Dxx=(gradient2(Dx,'x'));
Dxy=(gradient2(Dx,'y'));
clear Dx;

function D = gradient2(F,option)
% Example:
%
% Fx = gradient2(F,'x');

[k,l] = size(F);
D  = zeros(size(F),class(F));

switch lower(option)
case 'x'
    % Take forward differences on left and right edges
    D(1,:) = (F(2,:) - F(1,:));
    D(k,:) = (F(k,:) - F(k-1,:));
    % Take centered differences on interior points
    D(2:k-1,:) = (F(3:k,:)-F(1:k-2,:))/2;
case 'y'
    D(:,1) = (F(:,2) - F(:,1));
    D(:,l) = (F(:,l) - F(:,l-1));
    D(:,2:l-1) = (F(:,3:l)-F(:,1:l-2))/2;
otherwise
    disp('Unknown option')
end

function I=imgaussian(I,sigma,spacing,siz)
% IMGAUSSIAN filters an 1D, 2D color/greyscale or 3D image with an
% Gaussian filter. This function uses for filtering IMFILTER or if
% compiled the fast  mex code imgaussian.c . Instead of using a
% multidimensional gaussian kernel, it uses the fact that a Gaussian
% filter can be separated in 1D gaussian kernels.
%
% J=IMGAUSSIAN(I,SIGMA,SIZE)
%
% inputs,
%   I: 2D input image
%   SIGMA: The sigma used for the Gaussian kernel
%   SPACING: input image spacing
%   SIZ: Kernel size (single value) (default: sigma*6)
%
% outputs,
%   I: The gaussian filtered image
%

if(~exist('siz','var')), siz=sigma*6; end

if(sigma>0)

    % Filter each dimension with the 1D Gaussian kernels\
    x=-ceil(siz/spacing(1)/2):ceil(siz/spacing(1)/2);
    H = exp(-(x.^2/(2*(sigma/spacing(1))^2)));
    H = H/sum(H(:));
    Hx=reshape(H,[length(H) 1]);

    x=-ceil(siz/spacing(2)/2):ceil(siz/spacing(2)/2);
    H = exp(-(x.^2/(2*(sigma/spacing(2))^2)));
    H = H/sum(H(:));
    Hy=reshape(H,[1 length(H)]);

    I=imfilter(imfilter(I,Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate');
end

function [Lambda1,Lambda2,Ix,Iy]=eigvalOfHessian2D(Dxx,Dxy,Dyy)
% This function calculates the eigenvalues from the hessian matrix,
% sorted by abs value, and the eigenvector of Lambda1 (the
% smaller-magnitude eigenvalue = direction of minimal curvature = the
% direction running along a tubular structure). Adapted from D. Kroon's
% eig2image.m.

% Compute the eigenvectors of J, v1 and v2
tmp = sqrt((Dxx - Dyy).^2 + 4*Dxy.^2);
v2x = 2*Dxy; v2y = Dyy - Dxx + tmp;

% Normalize
mag = sqrt(v2x.^2 + v2y.^2); i = (mag ~= 0);
v2x(i) = v2x(i)./mag(i);
v2y(i) = v2y(i)./mag(i);

% The eigenvectors are orthogonal
v1x = -v2y;
v1y = v2x;

% Compute the eigenvalues
mu1 = 0.5*(Dxx + Dyy + tmp);
mu2 = 0.5*(Dxx + Dyy - tmp);

% Sort eigen values by absolute value abs(Lambda1)<abs(Lambda2)
check=abs(mu1)>abs(mu2);

Lambda1=mu1; Lambda1(check)=mu2(check);
Lambda2=mu2; Lambda2(check)=mu1(check);

% NOTE: v2 = (v2x, v2y) as constructed above is actually the eigenvector
% of mu1, and v1 = (v1x, v1y) (its orthogonal complement) is the
% eigenvector of mu2 -- so the pairing with Lambda1/Lambda2 below is
% swapped relative to the v1/v2 naming (verified numerically against a
% synthetic straight-line test image; a naive v1<->Lambda1 pairing comes
% out rotated 90 degrees).
Ix=v2x; Ix(check)=v1x(check);
Iy=v2y; Iy(check)=v1y(check);
