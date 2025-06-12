%% shaded errorbar
function ShadedErrorbar(meanData, errData,xCoor,Color)
% eg. Color.Shade='Blue'; Color.ShadeAlpha=0.5; Color.Line='Blue';
%%%%%%%% Example code %%%%%%%%%%%%%%%%%
% xCoor=[1:50];
% figure;
% Color.Shade='b';
% Color.ShadeAlpha=0.3;
% Color.Line='b';
% ShadedErrorbar(V',Vstd',xCoor,Color);
% Color.Shade='r';
% Color.ShadeAlpha=0.3;
% Color.Line='r';
% hold on
% ShadedErrorbar(VBB',VBBstd',xCoor,Color);
% xlabel('Vessel Index')
% ylabel('V [mm/s]')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on;
xVector=[xCoor,fliplr(xCoor)];
patch=fill(xVector,[meanData+errData, fliplr(meanData-errData)], Color.Shade);
set(patch,'edgecolor','none');
set(patch,'FaceAlpha',Color.ShadeAlpha);
hold on;
plot(xCoor, meanData,'color',Color.Line,'LineWidth',2);