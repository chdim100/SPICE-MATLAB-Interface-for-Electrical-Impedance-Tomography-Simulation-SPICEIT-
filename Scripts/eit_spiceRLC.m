function spice = eit_spiceRLC(img,path2lib,name,f,tissue_type,electrode_size)
%function spice = eit_spice(img, [name])
% Converts an EIT FEM model with assigned conductivities (an EIDORS "img") to a
% model reduced, fully connected mesh of resistors in SPICE format.
% If the FEM model has complex valued conductivities, the mesh will be an RLC
% mesh network.
%
% An optional subcircuit 'name' can be provided.
%
% DONE! complex value support (DIMAS 2020)
% TODO fix electrode ordering for mixed PEM/CEM electrodes
%
% CITATION_REQUEST:
% AUTHOR: A Boyle and A Adler
% TITLE: Integrating Circuit Simulation with EIT FEM Models
% JOURNAL: 19th International Conference on Biomedical Applications of Electrical Impedance Tomography, Edinburgh, UK
% YEAR: 2018
%

%  (C) 2018 A. Boyle, License: GPL version 2 or version 3

%%%%%%Inputs:
%%%-img: 3D FEM model
%%%-path2lib: path to SPICE library (where Netgen .cir models are stored)
%%%-name: 'name of the model'
%%%-f: input signal frequency (Hz)
%%%-tissue_type: 'Thorax' for thorax, 'Head' for head structure
%%%-electrode_size: '003' or '005' (normalized circular electrode radius)

if ischar(img) & strcmp(img, 'UNIT_TEST') unit_test(); end

if nargin == 1
    name = ['eit_',num2str(f),'Hz'];
end

Dprime = model_reduce(img);
%  disp(full(1./Dprime))
%  disp(full(-1./(Dprime - tril(Dprime))));
spice = netlist(Dprime,name,f);

if nargout == 0
        if strcmp(tissue_type,'Head')
            folder=[path2lib 'Structures\New_Structs_2020\' tissue_type '\'];
        else
            folder=[path2lib 'Structures\New_Structs_2020\' tissue_type '\' electrode_size '\'];
        end
    filename = [folder name '.cir' ];
    ind=1;
    while exist(filename,'file')==2
        ind=ind+1;
        filename = [folder name '_' num2str(ind) '.cir' ];
    end
    FILE = fopen(filename, 'w');
    fprintf(FILE,'%s\n',spice{:});
    fclose(FILE);
    eidors_msg(['saved SPICE netlist to ' filename]);
    return
end
end

function Dprime = model_reduce(img)
Y = calc_system_mat(img); Y= Y.E;
nn= num_nodes(img);
% Decompose into blocks; assumes that the nn+1:end nodes are CEM electrodes
rm = 1:nn; % nodes to fold
kp = nn+1:size(Y,1); % nodes to keep
% Now handle PEM electrodes, by transferring nodes between the rm and el sets
for i = 1:length(img.fwd_model.electrode)
    el = img.fwd_model.electrode(i);
    if length(el.nodes) == 1 %if electrodes are just points
        rm(rm == el.nodes) = []; %remove those points from rm
        kp(end+1) = el.nodes; %and add them to the electrode part
    end
end
% Note: C = B' ... we don't need to calculate it for symmetric matrices
A = Y(rm,rm); B= Y(rm,kp); D = Y(kp,kp);
%Dprime  = D - B'*inv(A)*B;
Dprime  = D - B'/A*B;
o=[];
end

function out = netlist(Dprime, name,f)
nn = size(Dprime,1);
ndr = floor(log10(nn*(nn-1)/2))+1; % number of digits for resistors
nde = floor(log10(nn))+1; % number of digits for electrodes
str = ['.subckt ' name ];
for ii = 1:nn
    str = [ str sprintf([' %' num2str(nde) 'd'], ii) ];
end
out = { str };
%str = ['R%' num2str(ndr) 'd %' num2str(nde) 'd %' num2str(nde) 'd %s'];
strR = ['R%d %' num2str(nde) 'd %' num2str(nde) 'd %s'];
strC = ['C%d %' num2str(nde) 'd %' num2str(nde) 'd %s'];
strL = ['L%d %' num2str(nde) 'd %' num2str(nde) 'd %s'];
rr=1;
cc=1;
ll=1;
for ii = 1:nn
    for jj = (ii+1):nn
        %         R=real(-1/Dprime(ii,jj));
        R=1/real(-Dprime(ii,jj));
        valR = sprintf('%3.12g',R);
        out(end+1,1) = { strrep(sprintf(strR,rr,ii,jj,valR),'+','') }; % we strip '+'
        rr = rr +1;
        %         CN=imag(-1/Dprime(ii,jj))/(2*pi*f);
        CN=1/(1/imag(-Dprime(ii,jj)));
        if CN>=0
            %             C=CN;
            C=CN/(2*pi*f);
            valC = sprintf('%3.12g',C);
            out(end+1,1) = { strrep(sprintf(strC,cc,ii,jj,valC),'+','') }; % we strip '+'
            cc=cc+1;
        else
            %             L=-CN;
            L=-1/(CN*2*pi*f);
            valL = sprintf('%3.12g',L);
            out(end+1,1) = { strrep(sprintf(strL,ll,ii,jj,valL),'+','') }; % we strip '+'
            ll=ll+1;
        end
    end
end
out(end+1,1) = { '.ends' };
end

function unit_test()
fmdl = mk_common_model('a2s',8);
img = mk_image(fmdl,1);
eit_spiceRC(img,'eit2')
end