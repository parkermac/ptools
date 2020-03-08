% makes a movie of local tide height

clear
close all

make_movie = 1;
if make_movie == 0
    ntt = 2;
else
    %ntt = 30;
    ntt = 25;
end

phome = '/Users/pm7/Documents/';
indir = 'LiveOcean_roms/output/cas6_v3_lo8b/f2019.07.04/';
infile = [phome,indir,'ocean_his_0001.nc'];
lon = nc_varget(infile,'lon_psi');
lat = nc_varget(infile,'lat_psi');
%%
figure;
set(gcf,'position',[100 100 1000.5 600]);

[NR, NC] = size(lon);
ttt = 0;
for tt = 1:ntt
    disp(['working on tt = ',num2str(tt)])
    tts = ['0000',num2str(tt)];
    infile = [phome,indir,'ocean_his_',tts(end-3:end),'.nc'];
    % plot ssh
    eta = squeeze(nc_varget(infile,'zeta'));
    etaa = NaN * ones(4, NR, NC);
    eata(1,:,:) = eta(1:end-1, 1:end-1);
    etaa(2,:,:) = eta(2:end, 1:end-1);
    etaa(3,:,:) = eta(2:end, 2:end);
    etaa(4,:,:) = eta(1:end-1, 2:end);
    etap = squeeze(nanmean(etaa, 1));
    
    h = surfl(lon,lat,etap);
    aa = axis;
    %axis([aa(1:4) -3 3])
    axis([-130 -122 42 52 -3 3])
    shading flat
    colormap copper
    lighting phong
    %set(gca,'DataAspectRatio',[1.46628 1 15]);
    set(gca,'DataAspectRatio',[1.46628 1 5]);
    %set(gca,'view',[-209 34]);
    set(gca,'view',[55 34]);
    set(gca,'fontsize',14);
    %set(gca,'view',[55 40]);
    
    %light('Position',[-125, -65, 10000],'Style','local')

    % add labels
    title('Sea Surface Height (m)','fontweight','bold','fontsize',20)
    %xlabel('Longitude (deg)')
    %ylabel('Latitude (deg)')
    
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