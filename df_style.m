function [canvas,colours,opt] = df_style(styID)
%% Description
%   Default plotting parameters specified by the style styID
%   Define additional styles within the switch statement. Users can create
%   their own styles for their projects and use the style package to point
%   towards the user-defined style
%
%
% Author
%   Naveed Ejaz (ejaz.naveed@gmail.com)

canvas  = 'blackonwhite';
opt     = [];
opt.save.journal        = 'brain';

switch(styID)
    case 'default'
        colours     = {'lightgreen','green','lightblue','blue'};
        opt.display.ax          = 'normal';
    case 'forcetrace'
        colours     = {'blue','green','red','magenta','black'};
        opt.general.markertype  = 'none';
        opt.display.ax          = 'square';
    case 'logslope'
        colours     = {'lightgray'};
        opt.general.markertype  = 'o';
        opt.general.markersize  = 4;
        opt.display.ax          = 'square';    
    case 'altenmuller_replication'
        colours     = {'orange','darkgray','lightgray'};
        opt.general.markertype  = 'o';
        opt.general.markersize  = 10;
        opt.legend.leglocation  = 'eastoutside';
        opt.display.ax          = 'square';            
    case 'group'
        colours     = {'blue','orange'};
        opt.general.markertype  = 'o';
        opt.general.markersize  = 10;
        opt.display.ax          = 'square';  
    case 'group_smallmarker'
        colours     = {'blue','orange'};
        opt.general.markertype  = 'o';
        opt.general.markersize  = 6;
        opt.display.ratio       = 'square';            
    case 'groupx3'
        colours     = {'blue','orange','yellow'};
        opt.general.markertype  = 'o';
        opt.general.markersize  = 6;
        opt.display.ax          = 'normal';            
end;

