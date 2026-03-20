function genSliderV2(Data)

fig=figure;
set(fig,'Toolbar','figure','NumberTitle','off')
% Create an axes to plot in
ax = axes('Position',[0.1300 0.1500 0.7750 0.8150]);
% sliders for epsilon and lambda
len = size(Data,3);
slider1_handle=uicontrol(fig,'Style','slider','Max',len,'Min',1,...
    'Value',1,'SliderStep',[1/(len-1) 10/(len-1)],...
    'Units','normalized','Position',[.02 .02 .9 .05]);
label = uicontrol(fig,'Style','text','Units','normalized','Position',[.02 .07 .14 .04],...
    'String','Choose frame');
% Set up callbacks
% climVals = [mean(min(Data,[],[1,2])),mean(max(Data,[],[1,2]))];
vars=struct('slider1_handle',slider1_handle,'Data',Data,'Axes',ax,'Label',label);%,'CLIMS',climVals);
set(slider1_handle,'Callback',{@slider1_callback,vars});
plotterfcn(vars)
% End of main file
end

% Callback subfunctions to support UI actions
function slider1_callback(~,~,vars)
    % Run slider1 which controls value of epsilon
    plotterfcn(vars)
end

function plotterfcn(vars)
    % Plots the image
    val = round(get(vars.slider1_handle,'Value'));
    imagesc(vars.Axes,vars.Data(:,:,val));
    colormap(vars.Axes,'gray');
    set(vars.Label,'String',num2str(val));
    axis equal; % added by Allen Zhou on 3/20/26
end