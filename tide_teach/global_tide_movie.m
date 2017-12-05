% makes a movie of global tide height

clear
close all


make_movie = 1;
if make_movie == 0
    ntt = 0;
else
    %ntt = 30;
    ntt = 359;
end

phome = '/Users/pm7/Documents/';
load([phome,'tools_data/forcing_data/tide/', ...
    'TPXO/Extractions/tpxo7p2_180.mat']);
x = Z.lon; y = Z.lat; z = Z.z;
[X,Y] = meshgrid(x,y);
mask = Z.mask; tmask = Z.topomask;

% plot selected constituents

fs1 = 16;
con = 'm2';
%con = 'k1';
con_str = [upper(con(1)),'_{',con(2),'}'];
eval(['amp = E_',con,'.amp;']);
eval(['phase = E_',con,'.phase;']);
amp(~mask) = NaN;
phase(~mask) = NaN; % 0-360

figure;
set(gcf,'position',[100 100 1000 600]);

ttt = 0;
for tt = 0:10:ntt
    disp(['working on tt = ',num2str(tt)])
    
    if 1
        surfl(x,y,amp .* cos(pi*(phase-tt)/180))
        aa = axis;
        axis([-180 180 -90 90 -3 3])
        shading flat
        colormap copper
        lighting phong
        set(gca,'DataAspectRatio',[1 1 .3]);
        set(gca,'view',[-11 32]);
        
        % add labels
        title([upper(con),' Height (m)'],'fontweight','bold')
        %xlabel('Longitude (deg)')
        %ylabel('Latitude (deg)')
    else
        pcolor(x,y,amp .* cos(pi*(phase-tt)/180));
        shading interp
        caxis([-.7, .7]);
        colorbar
        axis image
        title([con_str,' Tide Height (m) '],'fontsize',fs1+2,'fontweight','bold');
        xlabel('Longitude (deg)','fontsize',fs1)
        ylabel('Latitude (deg)','fontsize',fs1)
        set(gca,'fontsize',fs1);
        % add a coastline
        hold on
        contour(x,y,z,[-1, -1],'-k')
    end
    
    if make_movie % make a folder of jpegs for a movie
        if tt==0
            outdir = [phome,'ptools_output/tide/',con];
            if exist(outdir,'dir'); rmdir(outdir,'s'); end;
            mkdir(outdir);
        end
        set(gcf,'PaperPositionMode','auto');
        set(gcf,'position',[100 100 1000.5 600]);
        tts = ['0000',num2str(ttt)];
        print('-dpng',[outdir,'/plot_',tts(end-3:end),'.png']);
        ttt = ttt + 1;
        if tt<ntt; clf; end;
    end
    
end



