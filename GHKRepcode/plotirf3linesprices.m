
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make plots with 2 lines per panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot IRFs

% %% First figure
% 
% % titles for each panel
% titlelist = char('Output, CP (over a long horizon)','Output, CP (medium-run horizon)',...
%                  'Consumption, CP','Agg. Investment, CP',...
%                  'Consumption (output share)','Agg. Investment (output share)');
%              
% % labels for y axis in each panel
% ylabels = char('Percent Dev. From S.S.','Percent Dev. From S.S.','Percent Dev. From S.S.',...
%                'Percent Dev. From S.S.','P.Pt. Dev. From S.S.','P.Pt. Dev. From S.S.');
% 
% scale1 = 1;%/(f1_y_cp_irf(end)/f1_y_cp)/100;
% scale2 = 1;%/(f2_y_cp_irf(end)/f2_y_cp)/100;
% 
% % needs to match order in titlelist           
% second_line = 100*scale1*[f1_y_cp_irf/f1_y_cp, f1_y_cp_irf/f1_y_cp,...
%                   f1_c_cp_irf/f1_c_cp, f1_j_cp_irf/f1_j_cp,...
%                   f1_c_cp_share_irf,  f1_j_cp_share_irf ];%
%               
% third_line = 100*scale2*[f2_y_cp_irf/f2_y_cp, f2_y_cp_irf/f2_y_cp,...
%                    f2_c_cp_irf/f2_c_cp, f2_j_cp_irf/f2_j_cp,...
%                    f2_c_cp_share_irf,  f2_j_cp_share_irf ];
%                
% first_line = 100*scale1*[f3_y_cp_irf/f3_y_cp, f3_y_cp_irf/f3_y_cp,...
%                   f3_c_cp_irf/f3_c_cp, f3_j_cp_irf/f3_j_cp,...
%                   f3_c_cp_share_irf,  f3_j_cp_share_irf ];
% 
% % needs to match number of entries and order of variables included in
% % first_line (and second and third line if present).
% horizons = [nperiods_long,nperiods_short,...
%             nperiods_short,nperiods_short,...
%             nperiods_short,nperiods_short];
%         
% % plot figure
% makechartmix2(titlelist,legendlist,figtitle,first_line,second_line,third_line,horizons,ylabels)



%% Fifth figure, substitution wealth effect decomposition

deltaW = sum(beta.^(0:nperiods-1)'.*f1_c_irf/f1_c)*f1_c^(1-gamma);
deltaW1 = deltaW;
f1_wec_irf =  (1-beta)*deltaW/f1_c*ones(nperiods,1);
f1_sec_irf(1) = -beta/gamma*sum(beta.^(0:nperiods-1)'.*f1_r_irf);
for i = 2:nperiods
    f1_sec_irf(i) = f1_sec_irf(i-1)+1/gamma*f1_r_irf(i-1);
end
f1_sec_irf = f1_sec_irf';


deltaW = sum(beta.^(0:nperiods-1)'.*f2_c_irf/f2_c)*f2_c^(1-gamma);
deltaW2 = deltaW;
f2_wec_irf =  (1-beta)*deltaW/f2_c*ones(nperiods,1);
f2_sec_irf(1) = -beta/gamma*sum(beta.^(0:nperiods-1)'.*f2_r_irf);
for i = 2:nperiods
    f2_sec_irf(i) = f2_sec_irf(i-1)+1/gamma*f2_r_irf(i-1);
end
f2_sec_irf = f2_sec_irf';

deltaW = sum(beta.^(0:nperiods-1)'.*f3_c_irf/f3_c)*f3_c^(1-gamma);
deltaW3 = deltaW;
f3_wec_irf =  (1-beta)*deltaW/f3_c*ones(nperiods,1);
f3_sec_irf(1) = -beta/gamma*sum(beta.^(0:nperiods-1)'.*f3_r_irf);
for i = 2:nperiods
    f3_sec_irf(i) = f3_sec_irf(i-1)+1/gamma*f3_r_irf(i-1);
end
f3_sec_irf = f3_sec_irf';


% titles for each panel
titlelist = char('Output, CP (over a long horizon)','Output, CP (medium-run horizon)',...
                 'Consumption, CP','Agg. Investment, CP',...
                 'Wealth effect on Consumption','M Sector Output (share of aggregate)',...
                 'Substitution effect on Consumption','Relative Price of Equipment Investment');
             
% labels for y axis in each panel
ylabels = char('Percent Dev. From S.S.','Percent Dev. From S.S.',...
               'Percent Dev. From S.S.','Percent Dev. From S.S.',...
               'Percent Dev. From S.S.','Percent Dev. From S.S.',...
               'Percent Dev. From S.S.','Percent Dev. From S.S.');           
           
           
scale1 = 1;%/(f1_y_cp_irf(end)/f1_y_cp)/100;
scale2 = 1;%/(f2_y_cp_irf(end)/f2_y_cp)/100;
scale3 = 1;
% needs to match order in titlelist           
first_line = 100*scale1*[f1_y_cp_irf/f1_y_cp, f1_y_cp_irf/f1_y_cp,...
                  f1_c_cp_irf/f1_c_cp, f1_j_cp_irf/f1_j_cp,...
                  f1_wec_irf, f1_y_cp_share_h_irf,...
                  f1_sec_irf, f1_rpjepc_irf];
              
second_line = 100*scale2*[f2_y_cp_irf/f2_y_cp, f2_y_cp_irf/f2_y_cp,...
                   f2_c_cp_irf/f2_c_cp, f2_j_cp_irf/f2_j_cp,...
                   f2_wec_irf, f2_y_cp_share_h_irf,...
                   f2_sec_irf, f2_rpjepc_irf];
third_line = 100*scale3*[f3_y_cp_irf/f3_y_cp, f3_y_cp_irf/f3_y_cp,...
                   f3_c_cp_irf/f3_c_cp, f3_j_cp_irf/f3_j_cp,...
                   f3_wec_irf, f3_y_cp_share_h_irf,...
                   f3_sec_irf, f3_rpjepc_irf];
% needs to match number of entries and order of variables included in
% first_line (and second and third line if present).
horizons = [nperiods_long,nperiods_short,...
            nperiods_short,nperiods_short,...
            nperiods_short,nperiods_short,...
            nperiods_short,nperiods_long];


        
% plot figure
makechartmix2(titlelist,legendlist,figtitle,first_line,second_line,third_line,horizons,ylabels)

