% Nom des fichiers
nom_fic='escarp3.nc';
nom_fig='escarp3_hoev';

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


% CAS 1
% PENSER A COMMENTER LES BONNES LIGNES CI-DESSOUS POUR LE PLOT contourf !!!
% l'axe des X est l'axe des longitudes glamt
% celui des Y c'est le temps
% xphi0 c'est l'indice do la longitude ou on fait la coupe
xaxis=squeeze(glamt(1:end,1));
yaxis=squeeze(time(:));
xphi0 = 3;

% CAS 2
% PENSER A COMMENTER LES BONNES LIGNES CI-DESSOUS POUR LE PLOT contourf !!!
% l'axe des X est l'axe des longitudes glamt
% celui des Y c'est le temps
% yphi0 c'est l'indice de la latitude ou on fait la coupe
%xaxis=squeeze(time(:));
%yaxis=squeeze(gphit(1,1:end));
%yphi0 = 50;

figg = figure;
v=-0.02:0.002:0.02; % limites et espacements des contours

% CAS 1
hh=squeeze(hn(1:end,xphi0,1,1:end));
contourf(xaxis,yaxis,hh',v);
fig = [nom_fig '.xt']

% CAS 2
%hh=squeeze(hn(yphi0,1:end,1,1:end));
%contourf(xaxis,yaxis,hh,v);
%fig = [nom_fig '.yt']

colorbar;
shading flat;

% Choisir son format de sortie (tous les formats https://www.mathworks.com/help/matlab/ref/print.html)
print(figg,[fig '.jpg'],'-djpeg')
%print(figg,[fig '.png'],'-dpng')
%print(figg,[fig '.eps'],'-depsc')
%print(figg,[fig '.pdf'],'-dpdf')
