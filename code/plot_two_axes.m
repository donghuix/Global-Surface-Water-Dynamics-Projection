function [fig,axs,cb,pos2] = plot_two_axes(xv,yv,data1,data2,cmin,cmax,txt)
    
    cmap = flipud(blue2red(17));
    load('coastlines.mat');

    fig = figure; set(gcf,'Position',[10 10 1000 1200]);
    axs(1) = subplot(2,1,1);
    patch(xv,yv,data1,'LineStyle','none'); 
    clim([cmin cmax]);  hold on; grid on;
    ylim([-60 80]); xlim([-180 180]);%set(gca,'XTick',[],'YTick',[]); 
    plot(coastlon,coastlat,'k-','LineWidth',1);
    colormap(gca,cmap);
    
    axs(2) = subplot(2,1,2);
    patch(xv,yv,data2,'LineStyle','none'); 
    clim([cmin cmax]); hold on; grid on;
    ylim([-60 80]); xlim([-180 180]);%set(gca,'XTick',[],'YTick',[]); 
    plot(coastlon,coastlat,'k-','LineWidth',1);
    colormap(gca,cmap); cb = colorbar;
    
    axs(2).Position(2) = axs(2).Position(2) + 0.075;
    axs(1).Position(1) = axs(1).Position(1) - 0.1;
    axs(2).Position(1) = axs(2).Position(1) - 0.1;
    axs(1).Position(3) = axs(1).Position(3) - 0.1;
    axs(2).Position(3) = axs(2).Position(3) - 0.1;
    
    pos2 = get(axs(2),'Position');
    cb.Location     = 'south';
    cb.Position(1)  = pos2(1);
    cb.Position(2)  = pos2(2)-0.075;
    cb.Position(3)  = pos2(3);
    cb.Position(4)  = 0.02;
    cb.AxisLocation = 'out';
    cb.FontSize     = 13;
    title(cb,txt,'FontSize',15,'FontWeight','bold');
%     add_title(axs(1),'(a) SSP126',20,'out');
%     add_title(axs(2),'(b) SSP585',20,'out');
end