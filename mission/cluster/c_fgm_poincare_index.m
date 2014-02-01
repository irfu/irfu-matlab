function indices=c_fgm_poincare_index(vec1,vec2,vec3,vec4)
%
%Example:
%indices=c_fgm_poincare_index(B1,B2,B3,B4);
%B? is in any coordinates, 1st colomn is time, 2nd to 4th colomn is Bx,By,Bz
%
%------written by Yunhui Hu, Dec.20.2006 in Wuhan------------
%------modified by Shiyong Huang, May.17.2012 in IRF---------
n=size(vec1,2);
if n>3
    vec2=irf_resamp(vec2,vec1);
    vec3=irf_resamp(vec3,vec1);
    vec4=irf_resamp(vec4,vec1);
    time=vec1(:,1);
    vec1=vec1(:,2:4);
    vec2=vec2(:,2:4);
    vec3=vec3(:,2:4);
    vec4=vec4(:,2:4);
end
lx=size(vec1,1);
Map_sc1=zeros(lx,3);
Map_sc2=zeros(lx,3);
Map_sc3=zeros(lx,3);
Map_sc4=zeros(lx,3);
% project to 'M' space
   vl1= sqrt(dot(vec1,vec1,2));
   vl2= sqrt(dot(vec2,vec2,2));
   vl3= sqrt(dot(vec3,vec3,2));
   vl4= sqrt(dot(vec4,vec4,2));
   for i=1:3
    Map_sc1(:,i)=vec1(:,i)./vl1;
    Map_sc2(:,i)=vec2(:,i)./vl2;
    Map_sc3(:,i)=vec3(:,i)./vl3;
    Map_sc4(:,i)=vec4(:,i)./vl4;
    end
    %calculate the solid angle for each triangles
    area1=c_fgm_solidangle(Map_sc1,Map_sc2,Map_sc3);
    area2=c_fgm_solidangle(Map_sc1,Map_sc4,Map_sc2);
    area3=c_fgm_solidangle(Map_sc1,Map_sc3,Map_sc4);
    area4=c_fgm_solidangle(Map_sc2,Map_sc4,Map_sc3);
    indices=area1+area2+area3+area4;
%add time tags
if n>3
    indices=[time,indices];
else
    indices=indices;
end