function myfig(varargin)

mf.wwf                  =0.8;          % width-factor for plot window
mf.fig_fs               =20;           % [-] Figure Font Size
mf.whf                  =0.8;          % height-factor for plot window
mf.linw                 =5;            % LineWidth

mf.scrsz       =get(0,'ScreenSize');
mf.scrw        =mf.scrsz(3);
mf.scrh        =mf.scrsz(4);
mf.fwo         =(1-mf.wwf)*0.5;
mf.fho         =(1-mf.whf)*0.5;
mf.fwpos       =[mf.scrw*mf.fwo,mf.scrh*mf.fho,mf.scrw*mf.wwf,mf.scrh*mf.whf];

set(0,'DefaultLineLineWidth',2);
set(0,'DefaultaxesXGrid','on');
set(0,'DefaultaxesYGrid','on');
set(0,'DefaultaxesFontSize',mf.fig_fs);
set(0,'DefaultaxesFontName','Times New Roman');

if isempty(varargin)==1
    figure('Color','white','Position',mf.fwpos)
    hold on
else
    numfig=varargin{:};
    figure(numfig)
    a=get(numfig,'Color');
    if a(1)==1
    else
    set(gcf,'Color','white','Position',mf.fwpos)
    hold on
    grid on
    end
end

