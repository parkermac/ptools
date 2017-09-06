% makes a movie of local tide height

clear
close all

make_movie = 1;
if make_movie == 0
    ntt = 2;
else
    %ntt = 30;
    ntt = 73;
end

phome = '/Users/pm7/Documents/';
indir = 'LiveOcean_roms/output/cascadia1_base_lobio1/f2017.08.05/';
infile = [phome,indir,'ocean_his_0002.nc'];
lon = nc_varget(infile,'lon_rho');
lat = nc_varget(infile,'lat_rho');

figure;
set(gcf,'position',[100 100 1000.5 600]);

ttt = 0;
for tt = 2:ntt
    disp(['working on tt = ',num2str(tt)])
    tts = ['0000',num2str(tt)];
    infile = [phome,indir,'ocean_his_',tts(end-3:end),'.nc'];
    % plot ssh
    eta = squeeze(nc_varget(infile,'zeta'));
    surfl(lon,lat,eta);
    aa = axis;
    axis([aa(1:4) -3 3])
    %shading flat
    colormap copper
    lighting phong
    set(gca,'DataAspectRatio',[1.46628 1 15]);
    %set(gca,'view',[-209 34]);
    set(gca,'view',[55 34]);

    % add labels
    title('Sea Surface Height (m)','fontweight','bold')
    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')
    
    if make_movie % make a folder of jpegs for a movie
        if ttt==0
            outdir = [phome,'ptools_output/tide/LO'];
            if exist(outdir,'dir'); rmdir(outdir,'s'); end;
            mkdir(outdir);
        end
        set(gcf,'PaperPositionMode','auto');
        ttts = ['0000',num2str(ttt)];
        print('-dpng',[outdir,'/plot_',ttts(end-3:end),'.png']);
        ttt = ttt + 1;
        if tt<ntt; clf; end;
    end
    
end