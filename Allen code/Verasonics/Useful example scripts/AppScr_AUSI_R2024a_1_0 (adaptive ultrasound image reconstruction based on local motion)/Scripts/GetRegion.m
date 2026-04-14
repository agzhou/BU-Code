function GetRegion(imageBuff2)

% Copyright (C) 2001-2025, Verasonics, Inc.
% All worldwide rights and remedies under all intellectual
% property laws and industrial property laws are reserved.

   PData = evalin('base','PData');
   PdataREF = evalin('base','PdataREF');
   SE1 = evalin('base','SE1');
   dim1 = evalin('base','dim1');
   radius1 = evalin('base','radius1');
   fc = evalin('base','Trans.frequency');
   SoS = evalin('base','Resource.Parameters.speedOfSound');
   fPos = evalin('base','Resource.DisplayWindow.Position');

   mask = mean(squeeze(imageBuff2 ),3);
   mask = my_imclose(mask,SE1,dim1,radius1);
   mask = mask>median(mask);

   [pixels,~] = find(reshape(mask,[],1)~=0);
   PData(1).Region(end).PixelsLA =  int32(pixels);
   PData(1).Region(end).numPixels = length(PData(1).Region(end).PixelsLA);

   R2PixelsLA = PData(1).Region(end).PixelsLA;
   for jj = 2:length(PData(1).Region)-1
       newR1PixelsLA = PdataREF(1).Region(jj).PixelsLA;
       PData(1).Region(jj).PixelsLA  = setdiff(newR1PixelsLA,R2PixelsLA);
       PData(1).Region(jj).numPixels = length(PData(1).Region(jj).PixelsLA);
   end

   scaleToWvl = fc/(SoS/1000); 
   x = size(mask,2)*PData(1).PDelta(1)/ scaleToWvl; z = size(mask,1)*PData(1).PDelta(3)/ scaleToWvl;
   ff = figure(44);
   ff.Position = [(fPos(1)+fPos(3)+300) fPos(2)+50 fPos(3)-100 fPos(4)-100];  
   imagesc([-0.5*x 0.5*x],[0 z],mask),axis equal
   assignin('base','PData',PData);

   Control = evalin('base','Control');
   if isempty(Control(1).Command), n=1; else, n=length(Control)+1; end
   Control(n).Command = 'update&Run';
   Control(n).Parameters = {'PData','Recon'};
   assignin('base','Control', Control);
   assignin('base','mask',mask);
end

