  %--- Generate RF signals using a matrix transducer ---%
  %-- 3-MHz matrix array with 32x32 elements
  param = [];
  param.fc = 3e6;
  param.bandwidth = 70;
  param.width = 250e-6;
  param.height = 250e-6;
  %-- Position of the elements (pitch = 300 microns)
  pitch = 300e-6;
  [xe,ye] = meshgrid(((1:32)-16.5)*pitch);
  param.elements = [xe(:).'; ye(:).'];
  %-- Focus position
  x0 = 0; y0 = 0; z0 = 3e-2;
  %-- Transmit time delays using TXDELAY3
  dels = txdelay3(x0,y0,z0,param);
  %-- Create random scatterers
  N = 10000;
  x = 2*(rand(1,N)-0.5)*4e-3;
  y = 2*(rand(1,N)-0.5)*4e-3;
  z = rand(1,N)*6e-2;
  RC = hypot(rand(1,N),rand(1,N));
  %-- Simulate RF signals
  [RF,param] = simus3(x,y,z,RC,dels,param);
  %-- Display the elements and the scatterers
  figure
  scatter3(x*1e3,y*1e3,z*1e3,30,RC,'filled')
  colormap(cool)
  hold on
  scatter3(xe*1e3,ye*1e3,0*xe,3,'b','filled')
  axis equal, box on
  set(gca,'zdir','reverse')
  zlabel('[mm]')
  title([int2str(N) ' scatterers'])
  %-- Display the position of the elements
  figure, hold on
  plot(xe*1e3,ye*1e3,'ko','MarkerFaceColor','c')
  xlabel('x [mm]'), ylabel('y [mm]')
  axis square
  %-- Choose 4 elements randomly
  n = sort(randi(1024,1,4));
  plot(xe(n)*1e3,ye(n)*1e3,'ro','MarkerFaceColor','r')
  for k = 1:4
      text(xe(n(k))*1e3+0.3,ye(n(k))*1e3,int2str(n(k)),...
          'Color','r','BackgroundColor','w')
  end
  title('32{\times}32 elements')
  %-- Display their RF signals
  figure
  tl = tiledlayout(4,1);
  title(tl,'RF signals')
  RF = RF/max(RF(:,n),[],'all');
  for k = 1:4
      nexttile
      plot((0:size(RF,1)-1)/param.fs*1e6,RF(:,n(k)))
      title(['Element #' int2str(n(k))])
      ylim([-1 1])
  end