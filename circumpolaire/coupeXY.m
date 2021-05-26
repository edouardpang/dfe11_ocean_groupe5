% Nom des fichiers
nom_fic='test1.nc';
nom_fig='test1_XY';

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

% Choix des axes : ATTENTION, il faut enlever les bords des matrices
xaxis=squeeze(glamt(2:end-1,1));
yaxis=squeeze(gphit(1,2:end-1));

% choix du pas de temps où on fait la coupe XY
ntime0 = 3; 

figg = figure;
% limites et espacements des contours à ajuster en fonction des min max du champs à tracer
v=-0.02:0.002:0.02; 

% ATTENTION, il faut enlever les bords des matrices
hh=squeeze(hn(2:end-1,2:end-1,1,ntime0));
mmi = min(min(hh))
mma = max(max(hh))

contourf(xaxis,yaxis,hh',v);
fig = [nom_fig '.xy']

colorbar;
shading flat;

% Choisir son format de sortie (tous les formats https://www.mathworks.com/help/matlab/ref/print.html)
print(figg,[fig '.jpg'],'-djpeg')
%print(figg,[fig '.png'],'-dpng')
%print(figg,[fig '.eps'],'-depsc')
%print(figg,[fig '.pdf'],'-dpdf')
