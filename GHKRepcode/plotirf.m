
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot IRFs

%% First figure

% titles for each panel
titlelist = char('Output, CP','Output, CP',...
                 'Consumption, CP','Agg. Investment, CP',...
                 'Consumption (output share)','Agg. Investment (output share)');
             
% labels for y axis in each panel
ylabels = char('Percent Dev. From S.S.','Percent Dev. From S.S.','Percent Dev. From S.S.',...
               'Percent Dev. From S.S.','P.Pt. Dev. From S.S.','P.Pt. Dev. From S.S.');

scale0 = 1;%/(f0_y_cp_irf(end)/f0_y_cp)/100;

% needs to match order in titlelist           
first_line = 100*scale0*[f0_y_cp_irf/f0_y_cp, f0_y_cp_irf/f0_y_cp,...
                  f0_c_cp_irf/f0_c_cp, f0_j_cp_irf/f0_j_cp,...
                  f0_c_cp_share_irf,  f0_j_cp_share_irf ];
second_line = [];
third_line = [];

% needs to match number of entries and order of variables included in
% first_line (and second and third line if present).
horizons = [nperiods,nperiods_short,...
            nperiods_short,nperiods_short,...
            nperiods_short,nperiods_short];
        
% plot figure
makechartmix(titlelist,legendlist,figtitle,first_line,second_line,third_line,horizons,ylabels)


%% Second figure
    
% override the following line to change the horizon             
nnperiods = nperiods_short;

titlelist = char('L Structure Capital','H Structure Capital','L Equipment Capital','H Equipment Capital','L Labor','H Labor');
legendlist = char('');
figlabel = '';
makechart(titlelist,legendlist,figlabel,...
          100*scale0*[f0_k_sl_irf(1:nnperiods)/f0_k_sl, f0_k_sh_irf(1:nnperiods)/f0_k_sh,...
               f0_k_el_irf(1:nnperiods)/f0_k_el, f0_k_eh_irf(1:nnperiods)/f0_k_eh,...
               f0_n_l_irf(1:nnperiods)/f0_n_l,   f0_n_h_irf(1:nnperiods)/f0_n_h ])
      

%% Third figure
titlelist = char('L Structure Investment','H Structure Investment','L Equipment Investment','H Equipment Investment',...
    'Price of Equipment');
legendlist = char('');
figlabel = '';
makechart(titlelist,legendlist,figlabel,...
          100*scale0*[f0_j_sl_irf(1:nnperiods)/f0_j_sl, f0_j_sh_irf(1:nnperiods)/f0_j_sh,...
               f0_j_el_irf(1:nnperiods)/f0_j_el, f0_j_eh_irf(1:nnperiods)/f0_j_eh, f0_p_ce_irf(1:nnperiods)])
      

