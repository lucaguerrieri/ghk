
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make plots with 2 lines per panel
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

scale1 = 1;%/(f1_y_cp_irf(end)/f1_y_cp)/100;
scale2 = 1;%/(f2_y_cp_irf(end)/f2_y_cp)/100;

% needs to match order in titlelist           
first_line = 100*scale1*[f1_y_cp_irf/f1_y_cp, f1_y_cp_irf/f1_y_cp,...
                  f1_c_cp_irf/f1_y_cp, f1_j_cp_irf/f1_y_cp,...
                  f1_c_cp_share_irf,  f1_j_cp_share_irf ];%
              
second_line = 100*scale2*[f2_y_cp_irf/f2_y_cp, f2_y_cp_irf/f2_y_cp,...
                   f2_c_cp_irf/f2_y_cp, f2_j_cp_irf/f2_y_cp,...
                   f2_c_cp_share_irf,  f2_j_cp_share_irf ];
third_line = [];

% needs to match number of entries and order of variables included in
% first_line (and second and third line if present).
horizons = [nperiods_long,nperiods_short,...
            nperiods_short,nperiods_short,...
            nperiods_short,nperiods_short];
        
% plot figure
makechartmix(titlelist,legendlist,figtitle,first_line,second_line,third_line,horizons,ylabels)


%% Second figure
    
% override the following line to change the horizon             
nnperiods = nperiods_short;
titlelist = char('L Structure Capital','H Structure Capital','L Equipment Capital','H Equipment Capital','L Labor','H Labor',...
    'L Output','H Output');
figlabel = '';
line1 = 100*scale1*[f1_k_sl_irf(1:nnperiods)/f1_y_cp, f1_k_sh_irf(1:nnperiods)/f1_y_cp,...
               f1_k_el_irf(1:nnperiods)/f1_y_cp, f1_k_eh_irf(1:nnperiods)/f1_y_cp,...
               f1_n_l_irf(1:nnperiods)/f1_y_cp,   f1_n_h_irf(1:nnperiods)/f1_y_cp ,...
               f1_y_l_irf(1:nnperiods)/f1_y_cp,...
               f1_y_h_irf(1:nnperiods)/f1_y_cp...
               + f1_z_eh_irf(1:nnperiods)*f1_j_eh/f1_y_cp + f1_z_el_irf(1:nnperiods)*f1_j_el/f1_y_cp  ];

           
line2 = 100*scale2*[f2_k_sl_irf(1:nnperiods)/f2_y_cp, f2_k_sh_irf(1:nnperiods)/f2_y_cp,...
               f2_k_el_irf(1:nnperiods)/f2_y_cp, f2_k_eh_irf(1:nnperiods)/f2_y_cp,...
               f2_n_l_irf(1:nnperiods)/f2_y_cp,   f2_n_h_irf(1:nnperiods)/f2_y_cp,...
               f2_y_l_irf(1:nnperiods)/f2_y_cp,...
               f2_y_h_irf(1:nnperiods)/f2_y_cp...
               + f2_z_eh_irf(1:nnperiods)*f2_j_eh/f2_y_cp + f2_z_el_irf(1:nnperiods)*f2_j_el/f2_y_cp  ];

makechart(titlelist,legendlist,figlabel,line1, line2)
      

%% Third figure
titlelist = char('L Structure Investment','H Structure Investment','L Equipment Investment','H Equipment Investment',...
    'Price of Equipment');
figlabel = '';
line1 = 100*scale1*[f1_j_sl_irf(1:nnperiods)/f1_y_cp, f1_j_sh_irf(1:nnperiods)/f1_y_cp,...
               f1_j_el_irf(1:nnperiods)/f1_y_cp, f1_j_eh_irf(1:nnperiods)/f1_y_cp, f1_p_ce_irf(1:nnperiods)];
           
line2 = 100*scale2*[f2_j_sl_irf(1:nnperiods)/f2_y_cp, f2_j_sh_irf(1:nnperiods)/f2_y_cp,...
               f2_j_el_irf(1:nnperiods)/f2_y_cp + f2_z_el_irf(1:nnperiods)*f2_j_el/f2_y_cp,...
               f2_j_eh_irf(1:nnperiods)/f2_y_cp + f2_z_eh_irf(1:nnperiods)*f2_j_eh/f2_y_cp,...
               f2_p_ce_irf(1:nnperiods)];
           
           
makechart(titlelist,legendlist,figlabel,line1,line2)
      
%% Fourth figure, substitution wealth effect decomposition
titlelist = char('Wealth Effect on Consumption',...
                 'Substitution Effect on Consumption',...
                 'Overall Consumption Response')
figlabel = '';

deltaW = sum(beta.^(0:nperiods-1)'.*f1_c_irf/f1_c)*f1_c^(1-gamma);
f1_wec_irf =  (1-beta)*deltaW/f1_c*ones(nperiods,1);
f1_sec_irf(1) = -beta/gamma*sum(beta.^(0:nperiods-1)'.*f1_r_irf);
for i = 2:nperiods
    f1_sec_irf(i) = f1_sec_irf(i-1)+1/gamma*f1_r_irf(i-1);
end
f1_sec_irf = f1_sec_irf';

line1 = 100*scale1*[f1_wec_irf(1:nnperiods), f1_sec_irf(1:nnperiods),...
                    f1_c_irf(1:nnperiods)/f1_c];

deltaW = sum(beta.^(0:nperiods-1)'.*f2_c_irf/f2_c)*f2_c^(1-gamma);
f2_wec_irf =  (1-beta)*deltaW/f2_c*ones(nperiods,1);
f2_sec_irf(1) = -beta/gamma*sum(beta.^(0:nperiods-1)'.*f2_r_irf);
for i = 2:nperiods
    f2_sec_irf(i) = f2_sec_irf(i-1)+1/gamma*f2_r_irf(i-1);
end
f2_sec_irf = f2_sec_irf';
                
line2 = 100*scale2*[f2_wec_irf(1:nnperiods), f2_sec_irf(1:nnperiods),...
                    f2_c_irf(1:nnperiods)/f2_c];


makechart(titlelist,legendlist,figlabel,line1,line2)