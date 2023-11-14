clc
clear all
close all

foldernames = {
  '/Users/berends/Documents/Papers/UFEMISM2.0_part1/MISMIP/results_MISMIP_8km_spinup/'
  '/Users/berends/Documents/Papers/UFEMISM2.0_part1/MISMIP/results_MISMIP_8km_advance/'
  '/Users/berends/Documents/Papers/UFEMISM2.0_part1/MISMIP/results_MISMIP_8km_retreat/'
  };

wa = 600;
ha = 400;

margins_hor = [90,25];
margins_ver = [25,80];

H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);

xlabel( H.Ax{ 1,1},'Time (kyr)')
ylabel( H.Ax{ 1,1},'GL position (km)')

for fi = 1: length( foldernames)
  
  filename = [foldernames{ fi} '/main_output_ANT_grid.nc'];
  
  time = ncread( filename,'time');
  x_GL = zeros( size( time));
  x    = ncread( filename,'x');
  y    = ncread( filename,'y');
  
  imid = find( x==0);
  jmid = find( y==0);
  
  x = x( imid:end);
  
  for ti = 1: length( time)
    
    Hi = ncread( filename,'Hi',[imid,jmid,ti],[Inf,1,1]);
    Hb = ncread( filename,'Hb',[imid,jmid,ti],[Inf,1,1]);
    
    seawater_density = 1028;
    ice_density      = 910;
    TAF = Hi - max( 0.0, (0 - Hb) * (seawater_density / ice_density));
    
    for i = 1: length( x)-1
      if TAF( i) <= TAF( i+1)+1e-5
        TAF( 1:i) = TAF( 1:i) + 1e-3;
      end
    end
    
    x_GL( ti) = interp1( TAF,x,0);
    
  end
  
  line('parent',H.Ax{ 1,1},'xdata',time,'ydata',x_GL,'linewidth',3,'color','b');
end