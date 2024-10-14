% Need for Speed
% Caicedo-Pearce
% Untargeted moments from movers
% Fall 2024

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear, close all;
% clc, clear all;
%close all;

%Figure format
par=fig_format();

%Path
addpath('./mat/')
par.opt.dir_mat='./mat/';
par.opt.dir_fig='./figures/';
par.opt.dir_tab='./tables/';

%Saving options
par.opt.save_fig=0;
par.opt.save_tab=0;

%Running options
save_data=1;

warning('off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Upload data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Directory
dir='..\data\stata\'; 

%Strings to upload the data
data_str0='matched_pairs_1980_bottom_top_vs_bottom_';
data_opt_str='_c1_none_balanced_';
data_sf_str='gr1090_f';

%Cell with variables
data_varc={'lnnpat_f','speed','lnspeed','ln_speed_res','lf_cit3','lf_cit3_res','logxi','logxi_res'};
nvar=length(data_varc);

%Cell with regression specification
data_regc={'','wfe', 'strata','strata_wfe' };
nreg=length(data_regc);

for ivar=1:nvar
    for ireg=1:nreg
        %Upload data
        if strcmp(data_regc{ireg},'')
            filename =[dir data_str0 data_varc{ivar} data_opt_str data_regc{ireg} data_sf_str  '.csv'];
        else
            filename =[dir data_str0 data_varc{ivar} data_opt_str data_regc{ireg} '_' data_sf_str  '.csv'];
        end

        %Create data structure
        if strcmp(data_regc{ireg},'')
            data.(data_varc{ivar}).('base') = readtable(filename);
        else
            data.(data_varc{ivar}).(data_regc{ireg}) = readtable(filename);
        end
    end
end


%Save data
if save_data==1
        save([par.opt.dir_mat 'data_un_mom' data_opt_str data_sf_str '.mat'],'data')
end


%% No kmatch

%Directory
dir='..\data\stata\'; 

%Strings to upload the data
data_str0='move_no_kmatch_1980_bottom_top_vs_bottom_';
data_opt_str='_c1_none_balanced_';
data_sf_str='gr1090_f';

%Cell with variables
data_varc={'lnnpat_f','speed','lnspeed','ln_speed_res','lf_cit3','lf_cit3_res','logxi','logxi_res'};
nvar=length(data_varc);

%Cell with regression specification
data_regc={''};
nreg=length(data_regc);

for ivar=1:nvar
    for ireg=1:nreg
        %Upload data
        if strcmp(data_regc{ireg},'')
            filename =[dir data_str0 data_varc{ivar} data_opt_str data_regc{ireg} data_sf_str  '.csv'];
        else
            filename =[dir data_str0 data_varc{ivar} data_opt_str data_regc{ireg} '_' data_sf_str  '.csv'];
        end

        %Create data structure
        if strcmp(data_regc{ireg},'')
            data.(data_varc{ivar}).('base') = readtable(filename);
        else
            data.(data_varc{ivar}).(data_regc{ireg}) = readtable(filename);
        end
    end
end


%Save data
if save_data==1
        save([par.opt.dir_mat 'data_un_mom' data_str0 data_opt_str data_sf_str '.mat'],'data')
end