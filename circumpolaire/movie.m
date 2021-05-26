% Nom des fichiers
nom_fic='escarp3.nc';
nom_mov='escarp3';
nom_fic='standing0.001.low.nc';
nom_mov='standing0.001.low'; % le bon suffixe est mis par matlab
%VideoWriter.getProfiles   % pour avoir les codec disponibles
video = VideoWriter(nom_mov);
open(video);

% Lecture du fichier netcdf
time=ncread (nom_fic,'time'); 
hn=ncread (nom_fic,'HN');
un=ncread (nom_fic,'UN');
vn=ncread (nom_fic,'VN');
pv=ncread (nom_fic,'PV');
glamt=ncread (nom_fic,'GLAMT');
gphit=ncread (nom_fic,'GPHIT');
glamu=ncread (nom_fic,'GLAMU');
gphiu=ncread (nom_fic,'GPHIU');
glamv=ncread (nom_fic,'GLAMV');
gphiv=ncread (nom_fic,'GPHIV');

% creation des axes (le squeeze parce que glamt ... sont des tableaux 2D)
xaxis=squeeze(glamt(2:end-1,1));
yaxis=squeeze(gphit(1,2:end-1));

% creation de la figure
figure
v=0.0:0.002:0.012; % espacement contours pour contourf
%v=-0.1:0.05:0.6;
ntime=size(time)

% Boucle sur les pas de temps 
for j=1:5:ntime
%
  hh=squeeze(hn(2:end-1,2:end-1,1,j));
  hh=double(hh)';
  mmi = min(min(hh));
  mma = max(max(hh));
  %contourf(xaxis,yaxis,hh,v);
  surf(xaxis,yaxis,hh);view(0,90),shading flat;
  shading flat;
%
  %xtext=xaxis(1);
  %ytext=size(yaxis);ytext=yaxis(ytext(1)-1);
  %text(xtext,ytext+0.5,titre);
%
  caxis([-0.005,0.012]); % a modifier avec les bons min/max. Mieux de faire ca à la main
  %colorbar;
  axis equal; 
  axis([0 2000 0 2000]); % limite des axes x et y : [xmin xmax, ymin ymax]
%
  titre=['date: ' num2str(time(j)) ' (1/f units)'];
  %titre=['date: ' num2str(time(j)/86400.) ' days'];
  title(titre);
% creation video
  mov = getframe(1);
  writeVideo(video,mov);
%
end

close(video);
