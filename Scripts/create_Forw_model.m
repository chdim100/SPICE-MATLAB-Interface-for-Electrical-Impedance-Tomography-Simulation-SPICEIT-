function [img,elec_pos]=create_Forw_model(N, objects,symmetry,Rel,elec_pos)
%Rel=0.05;
if strcmp(symmetry,'Asymmetric')
    if nargin<5
    height=randn(1)*0.1+1;
    for i=1:length(objects)
        ox=objects{i}(:,1);
        oy=objects{i}(:,2);
        mx=mean(ox);
        my=mean(oy);
        oxm=ox-mx;
        oym=oy-my;
        oxmn=oxm*(0.1*randn(1)+1);
        oymn=oym*(0.1*randn(1)+1);
        oxmnm=oxmn+mx;
        oymnm=oymn+my;
        objectsnew{i}=[oxmnm oymnm];
        angles=90:-360/N:-270;
        neg=find(angles<0);
        angles(neg)=angles(neg)+360;
        angles=angles(1:end-1);
        angles_dev=angles.*(0.03*(rand(1,N)-0.5)+1);
        elec_z_dev=0.5*(0.05*randn(1,N)+1)+zeros(1,length(angles));
        elec_pos=[angles_dev; elec_z_dev]';
    end
    else
        height=1;
        objectsnew=objects;
        fprintf('electrode positions defined as input\n')
    end
    
else
    height=1;
    objectsnew=objects;
    elec_pos = [ N,                 % number of elecs per plane
                1,                   % equidistant spacing
                0.5]'; 
end
shape = { height,                      % height
          objectsnew, % contours
          %normally [4,50],
          [4,50],                 % perform smoothing with 50 points
          %1,
          0.04};                  % small maxh (fine mesh)
%normally: maxh==0.04
% elec_pos = [ N,                  % number of elecs per plane
%              1,                   % equidistant spacing
%              0.5]';               % a single z-plane


elec_shape = [Rel,               % radius
              0,                  % circular electrode
              0.01 ]';             % maxh (electrode refinement) 
%normally: maxh==0.01
fmdl = ng_mk_extruded_model(shape, elec_pos, elec_shape);
% this similar model is also available as:
% fmdl = mk_library_model('adult_male_16el_lungs');

[stim,meas_sel] = mk_stim_patterns(N,1,[0,1],[0,1],{'no_meas_current'}, 1);
fmdl.stimulation = stim;

img=mk_image(fmdl,1);
% img.elem_data(:)=0.35;
% img.elem_data(img.fwd_model.mat_idx{2})= 0.1812 + 0.0101i; % rlung
% img.elem_data(img.fwd_model.mat_idx{3})= 0.1812 + 0.0101i; % llung
% img.elem_data(img.fwd_model.mat_idx{4})= 0.1+0.002i; % spondilus
% img.elem_data(img.fwd_model.mat_idx{5})= 0.5+0.0025i; % heart

show_fem(img); view(0,70);
o=[];
end