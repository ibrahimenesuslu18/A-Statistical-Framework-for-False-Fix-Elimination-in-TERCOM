function [avg_score, best_score, best_pos, ratio,HATABUYUKflag, cands] = calculate_tercom_pos(pred_pos, terrain, P_INS, meas, varargin)
%TERCOM_MODULE Performs TERrain COntour Matching using Mean Absolute Difference.
% Inputs:
%   pred_pos        : [N x 2] predicted [lat, lon] values as degrees
%   P_INS           : [2 x 2] position unceratinty matrix of first prediction resolved in lat-lon
%   terrain         : terrain C interface object
%   errThreshold    : error threshold (m) for early termination
% Outputs:
%   bestPos         : [lat, lon] estimated start position in degrees
%   avgError        : Absolute error at best position match
%   score_map       : Score map for every test point

% Grid resolutions
skipSolution.MAD = 0;
skipSolution.MSD = 0;
skipSolution.DTW = 0;
skipSolution.PRO = 0;
resolution_m    = 30;
resolution_deg  = 1 / (3601);

meas = denoiseRadalt(meas);
% random_noise_bias = mvnrnd(zeros(1,2),P_INS,1);
random_noise_bias = mvnrnd(zeros(1,2),P_INS,1);

if ~isempty(varargin)
    H_first_point = varargin{1};

    H_bias = (H_first_point(1:2)  - pred_pos(1,1:2))+random_noise_bias;

    H_fp_relative = floor((H_first_point(1:2)  - pred_pos(1,1:2)-H_bias) / resolution_deg);  %#ok
    H_array = varargin{2};

    for i = 1:size(H_array,1)
        HoT_hybrid(i) = double(terrain.getMatrix(H_array(i,1),H_array(i,2),1));   %#ok
    end
end

% Number of trajectory points
N   = size(pred_pos, 1);
% Determine index offsets from first point
dRow = floor((pred_pos(:,1) - pred_pos(1,1)) / resolution_deg);
dCol = floor((pred_pos(:,2) - pred_pos(1,2)) / resolution_deg);


% Determine search window half-width from path span
spanRow = max(abs(dRow));
spanCol = max(abs(dCol));

%Determine covariance-based uncertainty window
sigma_lat = sqrt(P_INS(1, 1));
sigma_lon = sqrt(P_INS(2, 2));

% Convert to index space
idx_n = ceil(3 * sigma_lat / resolution_deg) + 2;
idx_e = ceil(3 * sigma_lon / resolution_deg) + 2;
% Determine search window size: cover entire trajectory

endPt = max([spanRow+idx_n, spanCol + idx_e]);
% matrixSize = ceil(2 * endPt/4) + 1;
matrixSize = ceil(2 * endPt) + 1;

% Center lat/lon for map
lat0 = pred_pos(1,1)+H_bias(1);
lon0 = pred_pos(1,2)+H_bias(2);

% lat0 = H_first_point(1);
% lon0 = H_first_point(2);

% Retrieve local DEM around predicted start (commented placeholder)
localDEM    = terrain.getMatrix(lat0, lon0, matrixSize);

tercomMap   = double(localDEM);

% tic;
% tercomMap = imresize(tercomMap,4,'bilinear'); 
% toc

% 
% [x_orig, y_orig] = meshgrid(1:size(tercomMap,1),1:size(tercomMap,2));
% newResCoef = 4;
% [xq, yq] = meshgrid(linspace(1,size(tercomMap,1),size(tercomMap,1)*newResCoef),linspace(1,size(tercomMap,1),size(tercomMap,1)*newResCoef));
% 
% % tic;
% tercomMap = interp2(x_orig,y_orig,tercomMap,xq,yq,'bilinear');
% toc

% figure()
% subplot(1,2,1)
% imagesc(Zq)
% subplot(1,2,2)
% imagesc(tercomMap)

% tercomMap   = flipud(tercomMap);

latspace    = lat0-endPt/3601 : 1/3601 : lat0+endPt/3601;
lonspace    = lon0-endPt/3601 : 1/3601 : lon0+endPt/3601;
mapSpace    = [latspace' lonspace'];

% Define the center index of the map
centerIdx   = ceil(size(tercomMap,2)/2); % initial guess index
centerIdx_scoreMap = [idx_n+1, idx_e+1];
% Construct the search points as latitude and longitude
searchArea  = [lat0 - idx_n/3601, lon0 - idx_e/3601;
               lat0 - idx_n/3601, lon0 + idx_e/3601;
               lat0 + idx_n/3601, lon0 + idx_e/3601;
               lat0 + idx_n/3601, lon0 - idx_e/3601];


% Core TERCOM search (nested loops with early exit)
score_map_MAD               = +inf(2*idx_n+1, 2*idx_e+1);
score_map_MSD               = +inf(2*idx_n+1, 2*idx_e+1);
score_map_TUKEY             = +inf(2*idx_n+1, 2*idx_e+1);
score_map_DTW               = +inf(2*idx_n+1, 2*idx_e+1);
% score_map_MI                = +inf(2*idx_n+1, 2*idx_e+1);
% score_map_ROC               = +inf(2*idx_n+1, 2*idx_e+1);
% score_map_CMB               = +inf(2*idx_n+1, 2*idx_e+1);


for i=-idx_n:1:idx_n
    for j=-idx_e:1:idx_e
        sumMAD = 0; sumMSD = 0; sumTukey=0; measHoT = zeros(1,N); mapHoT = zeros(1,N);
        for k=1:N
            rr = centerIdx + i + dRow(k);
            cc = centerIdx + j + dCol(k);

            sumMAD      = sumMAD + abs((pred_pos(k,3) - tercomMap(rr, cc)) - meas(k));
            sumMSD      = sumMSD + ((pred_pos(k,3) - tercomMap(rr, cc)) - meas(k))^2;

            measHoT(k)      = pred_pos(k,3) - meas(k);
            mapHoT(k)       = tercomMap(rr, cc);

            sumTukey = sumTukey +     calculateTukeyLoss(measHoT(k), mapHoT(k));
            % wi(k) = calculateTukeyLoss(measHoT(k), mapHoT(k));
        end
        % 
        % mapHoT_mean = mean(mapHoT);
        % measHoT_mean = mean(measHoT);

        % num = sum( (mapHoT-mapHoT_mean) .* (measHoT-measHoT_mean));
        % den = sqrt(sum( (mapHoT-mapHoT_mean).^2) * sum( (measHoT-measHoT_mean).^2)  );
        % 
        % % 
        % if den == 0
        %     sumNCC = 1;
        % else
        %     sumNCC = -(abs(num/den)-1);
        % end
        % 
        sumDTW = calculateDTWCost(measHoT, mapHoT);
        % 
        % sumCMB = log( ((sumDTW)^2 + (sumMAD)^2  + (sumMSD)) +1 );
        % sumTukey = wi * ((measHoT -mapHoT).^2)';
        
        
        
        score_map_MAD(idx_n + i +1, idx_e + j +1) = sumMAD/N;
        score_map_MSD(idx_n + i +1, idx_e + j +1) = sumMSD/N;
        score_map_TUKEY(idx_n + i +1, idx_e + j +1) = calculateMICost(measHoT, mapHoT)/N;
        score_map_DTW(idx_n + i +1, idx_e + j +1) = sumDTW/N;
        % score_map_MI(idx_n + i +1, idx_e + j +1)  = calculateMICost(measHoT, mapHoT)/N;
        % score_map_ROC(idx_n + i +1, idx_e + j +1) = diffHoT_cost(measHoT, mapHoT)/N;
        % score_map_CMB(idx_n + i +1, idx_e + j +1) = sumCMB/N;


    end

end




%%

% H_first_point_index =  H_fp_relative +centerIdx_scoreMap;


rNMS= 15;
cands_adaptive_MAD  = nms_botk_adaptive(score_map_MAD, rNMS);
cands_adaptive_MSD  = nms_botk_adaptive(score_map_MSD, rNMS);
cands_adaptive_NCC  = nms_botk_adaptive(score_map_TUKEY, rNMS);
cands_adaptive_DTW  = nms_botk_adaptive(score_map_DTW, rNMS);
% cands_adaptive_MI   = nms_botk_adaptive(score_map_MI , rNMS);
% cands_adaptive_ROC  = nms_botk_adaptive(score_map_ROC, rNMS);
% cands_adaptive_CMB  = nms_botk_adaptive(score_map_CMB, rNMS);

M_MAD = zeros(1,2);
M_MSD = zeros(1,2);
M_NCC = zeros(1,2);
M_DTW = zeros(1,2);
% M_MI  = zeros(1,2);
% M_ROC = zeros(1,2);
% M_CMB = zeros(1,2);

% totalCands = numel(cands_adaptive_MAD)/2;

for i = 1:2
    M_MAD(i) = score_map_MAD(cands_adaptive_MAD(i,1),cands_adaptive_MAD(i,2));
    M_MSD(i) = score_map_MSD(cands_adaptive_MSD(i,1),cands_adaptive_MSD(i,2));
    M_NCC(i) = score_map_TUKEY(cands_adaptive_NCC(i,1),cands_adaptive_NCC(i,2));
    M_DTW(i) = score_map_DTW(cands_adaptive_DTW(i,1),cands_adaptive_DTW(i,2));
    % M_MI (i) = score_map_MI(cands_adaptive_MI(i,1)  ,cands_adaptive_MI(i,2)) ;
    % M_ROC(i) = score_map_ROC(cands_adaptive_ROC(i,1),cands_adaptive_ROC(i,2));
    % M_CMB(i) = score_map_CMB(cands_adaptive_CMB(i,1),cands_adaptive_CMB(i,2));

end
%     neighScoresMAD = calculateCandsWeightedScore(score_map_MAD, cands_adaptive_MAD);
ratio.MAD = abs( (M_MAD(1)-M_MAD(2)) / M_MAD(2) )*100;
if ratio.MAD < 10
    skipSolution.MAD = 1;
end

ratio.MSD = abs( (M_MSD(1)-M_MSD(2)) / M_MSD(2) )*100;
if ratio.MSD < 10
    skipSolution.MSD = 1;
end
ratio.DTW = abs( (M_DTW(1)-M_DTW(2)) / M_DTW(2) )*100;
if ratio.DTW < 10
    skipSolution.DTW = 1;
end

ratio.MAD = abs( (M_MAD(1)-M_MAD(2)) / M_MAD(2) )*100;
ratio.MSD = abs( (M_MSD(1)-M_MSD(2)) / M_MSD(2) )*100;
ratio.NCC = abs( (M_NCC(1)-M_NCC(2)) / M_NCC(2) )*100;
ratio.DTW = abs( (M_DTW(1)-M_DTW(2)) / M_DTW(2) )*100;
% ratio.MI  = abs( (M_MI(1) -M_MI(2))  / M_MI(2)  )*100;
% ratio.ROC = abs( (M_ROC(1)-M_ROC(2)) / M_ROC(2) )*100;
% ratio.CMB = abs( (M_CMB(1)-M_CMB(2)) / M_CMB(2) )*100;


%     cands_adaptive_MAD(1,1),cands_adaptive_MAD(1,2)

% best position for mean absolute error
%     offset_abs_diff             = (bestLocal_abs_diff - [centerIdx, centerIdx]) * resolution_deg;
offset_MAD  = ([cands_adaptive_MAD(1,:)] - centerIdx_scoreMap(:)')*resolution_deg;
bestPos_MAD = [offset_MAD+H_bias 0]   + [pred_pos(end,:)];

offset_MSD  = ([cands_adaptive_MSD(1,:)] - centerIdx_scoreMap(:)')*resolution_deg;
bestPos_MSD = [offset_MSD 0]   + [pred_pos(end,:)];

offset_NCC  = ([cands_adaptive_NCC(1,:)] - centerIdx_scoreMap(:)')*resolution_deg;
bestPos_NCC = [offset_NCC 0]   + [pred_pos(end,:)];

offset_DTW  = ([cands_adaptive_DTW(1,:)] - centerIdx_scoreMap(:)')*resolution_deg;
bestPos_DTW = [offset_DTW+H_bias 0]   + [pred_pos(end,:)];


offset_MAD  = (cands_adaptive_MAD - centerIdx_scoreMap(:)')*resolution_deg;
cands.Pos_MAD = [offset_MAD, zeros(2,1)]   + [pred_pos(end,:)];


bestStartPos_DTW =  [offset_DTW+H_bias 0]   + [pred_pos(1,:)];


% 
% offset_MI  = ([cands_adaptive_MI(1,:)] - centerIdx_scoreMap(:)')*resolution_deg;
% bestPos_MI = [offset_MI+H_bias 0]   + [pred_pos(end,:)];
% 
% offset_ROC  = ([cands_adaptive_ROC(1,:)] - centerIdx_scoreMap(:)')*resolution_deg;
% bestPos_ROC = [offset_ROC+H_bias 0]   + [pred_pos(end,:)];
% 
% offset_CMB  = ([cands_adaptive_CMB(1,:)] - centerIdx_scoreMap(:)')*resolution_deg;
% bestPos_CMB = [offset_CMB+H_bias 0]   + [pred_pos(end,:)];

H_first_point_index =  H_fp_relative +centerIdx_scoreMap;

MAD_HATAidx = (cands_adaptive_MAD - H_first_point_index);
MAD_HATA = sqrt(MAD_HATAidx(:,1).^2+MAD_HATAidx(:,2).^2)*30;
DTW_HATAidx = (cands_adaptive_DTW - H_first_point_index);
DTW_HATA = sqrt(DTW_HATAidx(:,1).^2+DTW_HATAidx(:,2).^2)*30;
TKY_HATAidx = (cands_adaptive_NCC - H_first_point_index);
TKY_HATA = sqrt(TKY_HATAidx(:,1).^2+TKY_HATAidx(:,2).^2)*30;
  
roughness       = std(pred_pos(:,3)-meas);
roughnessDiff   = std(diff(pred_pos(:,3)-meas));

HATABUYUKflag.MAD = 0;
HATABUYUKflag.DTW = 0;
HATABUYUKflag.TKY = 0;

if min(MAD_HATA) > 120
    HATABUYUKflag.MAD = 1;
end
if min(DTW_HATA) > 120
    HATABUYUKflag.DTW = 1;
end
if min(TKY_HATA) > 120
    HATABUYUKflag.TKY = 1;
end

cands.MAD = cands_adaptive_MAD;
cands.DTW = cands_adaptive_DTW;
cands.TKY = cands_adaptive_NCC;
cands.merkez = centerIdx;
cands.Hibrit = H_first_point_index;
cands.lat0 = lat0 + (pred_pos(end,1)- pred_pos(1,1));
cands.lon0 = lon0 + (pred_pos(end,2)- pred_pos(1,2));
cands.scores = M_MAD';

%%


% figure();
% subplot(2,2,1)
% contourf(lonspace, latspace, tercomMap, 20); grid on;
% colormap('jet'); colorbar; hold on;
% xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
% title('Terrain Map and Vehicle Trajectories');
% % plot(lon_axis(round(true_traj_col)), lat_axis(round(true_traj_row)), 'g-', 'LineWidth', 2);
% plot([searchArea(:,2) ;searchArea(1,2)],[searchArea(:,1),;searchArea(1,1)],LineWidth=4,Color="black");
% scatter(H_first_point(2),H_first_point(1),180,"filled","pentagram",'r',"MarkerEdgeColor","black");
% plot(pred_pos(:,2),pred_pos(:,1), 'r:', 'LineWidth', 2);
% plot(H_array(:,2),H_array(:,1), 'black:', 'LineWidth', 2);
% 
% scatter(pred_pos(1,2),pred_pos(1,1),80,"filled","o",'r');
% scatter(bestPos_MAD(2),bestPos_MAD(1),80,"filled","square",'black');
% scatter(bestPos_DTW(2),bestPos_DTW(1),80,"filled","diamond",'black');
% 
% scatter(bestStartPos_DTW(2),bestStartPos_DTW(1),180,"filled","diamond","black")
% 
% subplot(2,2,3)
% imagesc(score_map_MAD);  xlim([1,2*idx_e+1]); ylim([1,2*idx_n+1]);hold on; grid on; colormap('jet'); colorbar; set(gca,'YDir','normal');
% scatter(cands_adaptive_MAD(1,2),cands_adaptive_MAD(1,1),80,"filled","square",'white');
% plot(cands_adaptive_MAD(:,2),cands_adaptive_MAD(:,1),'white+','MarkerSize',12,'LineWidth',1.5);
% scatter(H_first_point_index(2),H_first_point_index(1),180,"filled","pentagram",'r');
% 
% subplot(2,2,4)
% imagesc(score_map_DTW);  xlim([1,2*idx_e+1]); ylim([1,2*idx_n+1]);hold on; grid on; colormap('jet'); colorbar; set(gca,'YDir','normal');
% scatter(cands_adaptive_DTW(1,2),cands_adaptive_DTW(1,1),80,"filled","square",'white');
% plot(cands_adaptive_DTW(:,2),cands_adaptive_DTW(:,1),'white+','MarkerSize',12,'LineWidth',1.5);
% scatter(H_first_point_index(2),H_first_point_index(1),180,"filled","pentagram",'r');



%%


% plot(MAD_path_deg(:,2),MAD_path_deg(:,1), 'black:', 'LineWidth', 2);
% plot(MSD_path_deg(:,2),MSD_path_deg(:,1), 'black:', 'LineWidth', 2);
% legend('Terrain Elevation', 'True Position (GPS/INS)', 'INS-Only Trajectory', 'INS Start Position', 'MAD','MSD',"Search Region");




% dRowGorsel = ((pred_pos(:,1) - pred_pos(1,1)) / resolution_deg);
% dColGorsel = ((pred_pos(:,2) - pred_pos(1,2)) / resolution_deg);
% 
% MAD_path_index  = cands_adaptive_MAD(1,:)    + [(dRowGorsel) , (dColGorsel)];
% MAD_path_index = [MAD_path_index,800 * ones(N,1)];
% 
% MAD_path_index  = H_first_point_index(1,:)    + [(dRowGorsel) , (dColGorsel)];
% MAD_path_index = [MAD_path_index,800 * ones(N,1)];
% 
% vvv = size(tercomMap,1);
% xv = 1:vvv;
% yv = 1:vvv;
% [X,Y] = meshgrid(xv,yv);
% 
% 
% 
% figure()
% surf(X,Y,tercomMap,'EdgeColor','none');  hold on; grid on; colorbar; set(gca,'YDir','normal'); title("Map");
% scatter3(H_first_point_index(2),H_first_point_index(1),800,180,"filled","pentagram",'r');
% plot3(MAD_path_index(:,2),MAD_path_index(:,1),MAD_path_index(:,3),'black','LineWidth',3);
% 
% 
% for i = 1:size(MAD_path_index,1)
%     % x = MAD_path_index(i)
%     plot3([MAD_path_index(i,2) MAD_path_index(i,2)], [MAD_path_index(i,1) MAD_path_index(i,1)],[MAD_path_index(i,3), 0], 'r-','LineWidth',2)
% end
% % zlim([0 1200])
% % xlim([0 25])
% % ylim([0 25])
% 
% mask = false(size(tercomMap,1)-1,size(tercomMap,2)-1 );
% 
% GORSEL_index =  MAD_path_index + [2*ones(N,1),3*ones(N,1), zeros(N,1) ]; 
% GORSEL_index(11,:) = [6,11,800]; 
% GORSEL_index(14,:) = [6,11,800]; 
% GORSEL_index(15,:) = [6,11,800]; 
% GORSEL_index(21,:) = [6,11,800]; 
% GORSEL_index(24,:) = [6,11,800]; 
% GORSEL_index(25,:) = [6,11,800]; 
% GORSEL_index(27,:) = [6,11,800]; 
% GORSEL_index(28,:) = [6,11,800]; 
% GORSEL_index(34,:) = [6,11,800]; 
% GORSEL_index(35,:) = [6,11,800]; 
% GORSEL_index(37,:) = [6,11,800]; 
% GORSEL_index(38,:) = [6,11,800]; 
% GORSEL_index(31,:) = [6,11,800]; 
% 
% 
% for i = 1:size(GORSEL_index,1)
%     ix = find(xv <= GORSEL_index(i,1),1,'last');
%     iy = find(yv <= GORSEL_index(i,2),1,'last');
% 
%     if ~isempty(ix) && ix < numel(xv) && ~isempty(iy) && iy < numel(yv)
%         mask(iy,ix) = true;
%     end
% end
% surface('XData', X, 'YData', Y, 'ZData', tercomMap+1e-3, ...
%     'FaceColor',[1 0 0], ...
%     'FaceAlpha', 'flat', ...
%     'AlphaData', 0.65*mask);
% 
% 

% 
% if MAD_HATA > 90 || DTW_HATA > 90 || TKY_HATA > 90
% % if TKY_HATA > 30 && roughness > 0.0001
    % figure()
    % % subplot(2,3,1)
    % imagesc(score_map_MAD);  xlim([1,2*idx_e+1]); ylim([1,2*idx_n+1]);hold on; grid on; colormap('jet'); colorbar; set(gca,'YDir','normal');
    % scatter(cands_adaptive_MAD(1,2),cands_adaptive_MAD(1,1),80,"filled","square",'white');
    % plot(cands_adaptive_MAD(:,2),cands_adaptive_MAD(:,1),'white+','MarkerSize',12,'LineWidth',1.5);
    % scatter(H_first_point_index(2),H_first_point_index(1),180,"filled","pentagram",'r');
    % scatter(centerIdx_scoreMap(2),centerIdx_scoreMap(1),109080,"o",'white');
    % scatter(centerIdx_scoreMap(2),centerIdx_scoreMap(1),209080,"o",'white');
    % scatter(centerIdx_scoreMap(2),centerIdx_scoreMap(1),309080,"o",'white');
    % 
    % title("MAD score map");
    % 

    % subplot(2,3,2)
    % imagesc(score_map_DTW);  xlim([1,2*idx_e+1]); ylim([1,2*idx_n+1]);hold on; grid on; colormap('jet'); colorbar; set(gca,'YDir','normal');
    % scatter(cands_adaptive_DTW(1,2),cands_adaptive_DTW(1,1),80,"filled","square",'white');
    % plot(cands_adaptive_DTW(:,2),cands_adaptive_DTW(:,1),'white+','MarkerSize',12,'LineWidth',1.5);
    % scatter(H_first_point_index(2),H_first_point_index(1),180,"filled","pentagram",'r');
    % title("DTW score map");
    % 
    % subplot(2,3,3)
    % imagesc(score_map_TUKEY);  xlim([1,2*idx_e+1]); ylim([1,2*idx_n+1]);hold on; grid on; colormap('jet'); colorbar; set(gca,'YDir','normal');
    % scatter(cands_adaptive_NCC(1,2),cands_adaptive_NCC(1,1),80,"filled","square",'white');
    % plot(cands_adaptive_NCC(:,2),cands_adaptive_NCC(:,1),'white+','MarkerSize',12,'LineWidth',1.5);
    % scatter(H_first_point_index(2),H_first_point_index(1),180,"filled","pentagram",'r');
    % title("TUKEY score map");
    % 
    % subplot(2,3,5)
    % imagesc(tercomMap);  hold on; grid on; colormap('jet'); colorbar; set(gca,'YDir','normal');
    % title("Map");
% 
%     % mp_cnv = (size(tercomMap) - size(score_map_MAD))/2; 
% 
%     MAD_path_index  = (cands_adaptive_MAD(1,:))   + [(dRow) , (dCol)];
%     for sfjkld = 1: size(MAD_path_index,1)
%         % - centerIdx_scoreMap(:)'
%         HoT_harita(sfjkld) = tercomMap(MAD_path_index(sfjkld,1)- centerIdx_scoreMap(1)+centerIdx ,MAD_path_index(sfjkld,2)- centerIdx_scoreMap(2) +centerIdx);
%     end
%     subplot(2,3,6)
%     plot(pred_pos(:,3) - meas(:),'LineWidth',2);hold on; grid on;
%     plot(HoT_harita,'LineWidth',2);hold on; grid on;
%     legend("Real HoT","DTED HoT");
%     fprintf("Hata          = %.3f\nRoughness     = %.3f\nRoughnessDiff = %.3f\nRatioMAD      = %.3f\nRatioDTW      = %.3f\nRatioTUKEY    = %.3f\n",TKY_HATA,roughness,roughnessDiff,ratio.MAD,ratio.DTW,ratio.NCC);
%     close all;
% 
% end



best_pos = struct();
best_pos.MAD  = bestPos_MAD;
best_pos.MSD  = bestPos_MSD;
best_pos.NCC  = bestPos_NCC;
best_pos.DTW  = bestPos_DTW;
% best_pos.MI   = bestPos_MI;
% best_pos.ROC  = bestPos_ROC;
% best_pos.CMB  = bestPos_CMB;

avg_score = struct();
avg_score.MAD  = M_MAD(1);
avg_score.MSD  = M_MSD(1);
avg_score.NCC  = M_NCC(1);
avg_score.DTW  = M_DTW(1);
% avg_score.MI   = M_MI (1);
% avg_score.ROC  = M_ROC(1);
% avg_score.CMB  = M_CMB(1);


best_score = struct();
best_score.MAD  = M_MAD(1)/N;
best_score.MSD  = M_MSD(1)/N;
best_score.NCC  = M_NCC(1)/N;
best_score.DTW  = M_DTW(1)/N;
% best_score.MI   = M_MI (1)/N;
% best_score.ROC  = M_ROC(1)/N;
% best_score.CMB  = M_CMB(1)/N;


% %     DTW_dom = ratio.DTW > 2*ratio.MAD;
% %     MAD_dom = ratio.MAD > 2*ratio.DTW;
% 
% %     if DTW_dom
% %         best_pos.PRO = bestPos_DTW;
% %     elseif MAD_dom
% %         best_pos.PRO = bestPos_MAD;
% %     else
% %         best_pos.PRO = bestPos_DTW;
% %     end
% if (ratio.MAD >= ratio.DTW)
%     ratio.PRO = ratio.MAD;
%     best_pos.PRO = bestPos_MAD;
%     skipSolution.PRO = skipSolution.MAD;
% else
%     ratio.PRO = ratio.DTW;
%     best_pos.PRO = bestPos_DTW;
%     skipSolution.PRO = skipSolution.DTW;
%     %     else
%     %         ratio.PRO = ratio.MSD;
%     %         best_pos.PRO = bestPos_MSD;
%     %         skipSolution.PRO = skipSolution.MSD;
% end


%% outputs



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% map features
% 
% % TRI
% % terrainRuggednessIndex = sum(sum(abs(diff(tercomMap,1,1)))) + sum(sum(abs(diff(tercomMap,1,2))));
% [TRImean, TRIARSclassIndex] = compute_ARS_from_MAP(tercomMap);
% 
% 
% % TPI
% cellSize = 30;
% innerSize = min((2*idx_n +1),(2*idx_e +1));
% scSmall = 300;
% scLarge = 2100;
% [tpiLE, tpiPVS, tpiSmax, tpiSAR] = tpi_feature_pack(tercomMap, cellSize, innerSize, scSmall, scLarge);
% 
% win = 5;
% localMean = conv2(tercomMap, ones(win)/win^2 , 'same');
% topographicPositionIndex = mean(tercomMap(:) - localMean(:));
% 
% % sloce Aspect 2d gradient
% [dzdx, dzdy] = gradient(tercomMap,resolution_m);
% slope = atan(sqrt(dzdx.^2 + dzdy.^2));
% % aspect = atan2(dzdy, dzdx);
% meanSlopeMap = mean(slope(:));
% 
% % Plan & Profile Curvature
% [dxx, dxy] = gradient(dzdx);
% [~, dyy] = gradient(dzdy);
% planCurv = dxx + dyy;
% profileCurv = (dxx .* dzdx.^2 + 2*dxy.*dzdx.*dzdy + dyy.*dzdy.^2) ...
%     ./ (dzdx.^2 + dzdy.^2 + 1e-6);
% meanPlanCurvMap = mean(planCurv(:));
% meanProfileCurvMap = mean(profileCurv(:));
% 
% % histogram momentleri
% skewnessMap = skewness(tercomMap,0,'all');
% kurtosisMap = kurtosis(tercomMap,0,'all');
% 
% % yerel doruk says
% [nPeaksMap, nPeakDensityMap]= count_peaks(tercomMap);
% 
% % entropi
% % entropyMap = entropy(uint8(mat2gray(tercomMap)*255));
% entropyMap = my_entropy(tercomMap,255);
% %spread, variance
% varianceMap = var(tercomMap(:));
% 
% %
% meanMap                = mean(tercomMap,'all');
% minMap = min(tercomMap(:));
% maxMap = max(tercomMap(:));
% rangeMap =maxMap - minMap;
% iqrMap = iqr(tercomMap(:));
% 
% % 2D FFT
% F = fft2(tercomMap);
% Fs = fftshift(F);
% magFs = abs(Fs);
% % totalEnergy = sum(magF(:).^2);
% powFs = (abs(Fs).^2) / numel(tercomMap);
% FFTMaxPowFsMap = max(powFs,[],'all');
% FFTMeanPowFsMap = mean(powFs,'all');
% FFTMedianPowFsMap = median(powFs,'all');
% FFTstdPowFsMap = std(powFs,[],'all');
% 
% [M,N] = size(tercomMap);
% [u,v] = meshgrid(1:N,1:M);
% centroid_u = sum(sum(u.*magFs)) / sum(magFs(:));
% centroid_v = sum(sum(v.*magFs)) / sum(magFs(:));
% spectralCentroid = [centroid_v, centroid_u];
% 
% bw_u = sqrt(sum(sum(((u-centroid_u).^2).*magFs)) / sum(magFs(:)) );
% 
% FFTSpectralEntrpoyMap = spectral_entropy_radial(tercomMap);
% 
% 
% mapFeatures = struct();
% mapFeatures.mapTRImean                     = TRImean;
% mapFeatures.mapTRIARSclassIndex            = TRIARSclassIndex;
% mapFeatures.mapTPImean                     = topographicPositionIndex;
% mapFeatures.mapTPILE                       = tpiLE;
% mapFeatures.mapTPIPVS                      = tpiPVS;
% mapFeatures.mapTPISmax                     = tpiSmax;
% mapFeatures.mapTPISAR                      = tpiSAR;
% mapFeatures.mapMeanSlope                   = meanSlopeMap;
% mapFeatures.mapMeanPlanaryCurvature        = meanPlanCurvMap;
% mapFeatures.mapMeanProfileCurvature        = meanProfileCurvMap;
% mapFeatures.mapSkewness                    = skewnessMap;
% mapFeatures.mapKurtosis                    = kurtosisMap;
% mapFeatures.mapNPeaks                      = nPeaksMap;
% mapFeatures.mapNPeakDensity                = nPeakDensityMap;
% mapFeatures.mapEntropy                     = entropyMap;
% mapFeatures.mapMean                        = meanMap;
% mapFeatures.mapVariance                    = varianceMap;
% mapFeatures.mapMin                         = minMap;
% mapFeatures.mapMax                         = maxMap ;
% mapFeatures.mapRange                       = rangeMap;
% mapFeatures.mapIqr                         = iqrMap;
% mapFeatures.mapFFTspectralCentroid1        = spectralCentroid(1);
% mapFeatures.mapFFTspectralCentroid2        = spectralCentroid(2);
% mapFeatures.mapFFTbw_u                     = bw_u;
% mapFeatures.mapFFTRadialSpectralEntropy    = FFTSpectralEntrpoyMap;
% mapFeatures.mapFFTMaxPowFs                 = FFTMaxPowFsMap;
% mapFeatures.mapFFTMeanPowFs                = FFTMeanPowFsMap;
% mapFeatures.mapFFTMedianPowFs              = FFTMedianPowFsMap;
% mapFeatures.mapFFTstdPowFs                 = FFTstdPowFsMap;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% score map features
% scoreMapSize = size(score_map_MAD,2);
% centerScoreMap = floor(scoreMapSize/2);
% patchSizeScoreMap = 40;
% halfPatchSize = floor(patchSizeScoreMap/2);
% 
% % S_patch = score_map_MAD((centerScoreMap - halfPatchSize):(centerScoreMap + halfPatchSize+1),(centerScoreMap - halfPatchSize):(centerScoreMap + halfPatchSize+1));
% 
% meanScoreMap                = mean(score_map_MAD,'all');
% skewnessScoreMap            = skewness(score_map_MAD,0,'all');
% kurtosisScoreMap            = kurtosis(score_map_MAD,0,'all');
% % entropyScoreMap             = entropy(uint8(mat2gray(score_map_MAD)*255));
% entropyScoreMap = my_entropy(score_map_MAD,255);
% 
% [nPeaksScoreMap, nPeakDensityScoreMap]= count_peaks(score_map_MAD);
% 
% % sloce Aspect 2d gradient
% [dzdx, dzdy] = gradient(score_map_MAD,resolution_m);
% slope = atan(sqrt(dzdx.^2 + dzdy.^2));
% % aspect = atan2(dzdy, dzdx);
% meanSlopeScoreMap = mean(slope(:));
% 
% % Plan & Profile Curvature
% [dxx, dxy] = gradient(dzdx);
% [~, dyy] = gradient(dzdy);
% planCurvScoreMap = dxx + dyy;
% profileCurvScoreMap = (dxx .* dzdx.^2 + 2*dxy.*dzdx.*dzdy + dyy.*dzdy.^2) ...
%     ./ (dzdx.^2 + dzdy.^2 + 1e-6);
% meanPlanCurvScoreMap = mean(planCurvScoreMap(:));
% meanProfileCurvScoreMap = mean(profileCurvScoreMap(:));
% 
% % 2D FFT
% F = fft2(score_map_MAD);
% Fs = fftshift(F);
% magFs = abs(Fs);
% 
% powFs = (abs(Fs).^2) / numel(score_map_MAD);
% FFTMaxPowFsScoreMap = max(powFs,[],'all');
% FFTMeanPowFsScoreMap = mean(powFs,'all');
% FFTMedianPowFsScoreMap = median(powFs,'all');
% FFTstdPowFsScoreMap = std(powFs,[],'all');
% 
% [M,N] = size(score_map_MAD);
% [u,v] = meshgrid(1:N,1:M);
% centroid_u = sum(sum(u.*magFs)) / sum(magFs(:));
% centroid_v = sum(sum(v.*magFs)) / sum(magFs(:));
% spectralCentroid = [centroid_v, centroid_u];
% 
% bw_u = sqrt(sum(sum(((u-centroid_u).^2).*magFs)) / sum(magFs(:)) );
% 
% FFTSpectralEntrpoyScoreMap = spectral_entropy_radial(score_map_MAD);
% 
% varianceScoreMap = var(score_map_MAD(:));
% 
% %
% minScoreMap = min(score_map_MAD(:));
% maxScoreMap = max(score_map_MAD(:));
% rangeScoreMap =maxScoreMap - minScoreMap;
% iqrScoreMap = iqr(score_map_MAD(:));
% 
% 
% scoreMapFeatures = struct();
% 
% scoreMapFeatures.scoreMapMeanSlope                   = meanSlopeScoreMap;
% scoreMapFeatures.scoreMapMeanPlanaryCurvature        = meanPlanCurvScoreMap;
% scoreMapFeatures.scoreMapMeanProfileCurvature        = meanProfileCurvScoreMap;
% 
% scoreMapFeatures.scoreMapSkewness                    = skewnessScoreMap;
% scoreMapFeatures.scoreMapKurtosis                    = kurtosisScoreMap;
% scoreMapFeatures.scoreMapEntropy                     = entropyScoreMap;
% scoreMapFeatures.scoreMapNPeaks                      = nPeaksScoreMap;
% scoreMapFeatures.scoreMapNPeakDensity                = nPeakDensityScoreMap;
% 
% scoreMapFeatures.scoreMapMean                        = meanScoreMap;
% scoreMapFeatures.scoreMapVariance                    = varianceScoreMap;
% scoreMapFeatures.scoreMapMin                         = minScoreMap;
% scoreMapFeatures.scoreMapMax                         = maxScoreMap;
% scoreMapFeatures.scoreMapRange                       = rangeScoreMap;
% scoreMapFeatures.scoreMapIqr                         = iqrScoreMap;
% 
% scoreMapFeatures.scoreMapFFTspectralCentroid1        = spectralCentroid(1);
% scoreMapFeatures.scoreMapFFTspectralCentroid2        = spectralCentroid(2);
% scoreMapFeatures.scoreMapFFTbw_u                     = bw_u;
% scoreMapFeatures.scoreMapFFTSpectralEntrpoy          = FFTSpectralEntrpoyScoreMap;
% scoreMapFeatures.scoreMapFFTMaxPowFs                 = FFTMaxPowFsScoreMap;
% scoreMapFeatures.scoreMapFFTMeanPowFs                = FFTMeanPowFsScoreMap;
% scoreMapFeatures.scoreMapFFTMedianPowFs              = FFTMedianPowFsScoreMap;
% scoreMapFeatures.scoreMapFFTstdPowFs                 = FFTstdPowFsScoreMap;


    % H_first_point_index =  H_fp_relative +centerIdx_scoreMap;




% % % 
% figure()
% subplot(2,2,1)
% imagesc(score_map_MAD);  xlim([1,2*idx_e+1]); ylim([1,2*idx_n+1]);hold on; grid on; colormap('jet'); colorbar; set(gca,'YDir','normal');
% scatter(cands_adaptive_MAD(1,2),cands_adaptive_MAD(1,1),80,"filled","square",'white');
% plot(cands_adaptive_MAD(:,2),cands_adaptive_MAD(:,1),'white+','MarkerSize',12,'LineWidth',1.5);
% scatter(H_first_point_index(2),H_first_point_index(1),180,"filled","pentagram",'r');
% title("MAD score map");
% 
% % 
% % subplot(3,3,2)
% % imagesc(score_map_MSD);  xlim([1,2*idx_e+1]); ylim([1,2*idx_n+1]);hold on; grid on; colormap('jet'); colorbar; set(gca,'YDir','normal');
% % scatter(cands_adaptive_MSD(1,2),cands_adaptive_MSD(1,1),80,"filled","square",'white');
% % plot(cands_adaptive_MSD(:,2),cands_adaptive_MSD(:,1),'white+','MarkerSize',12,'LineWidth',1.5);
% % 
% subplot(2,2,2)
% imagesc(score_map_TUKEY);  xlim([1,2*idx_e+1]); ylim([1,2*idx_n+1]);hold on; grid on; colormap('jet'); colorbar; set(gca,'YDir','normal');
% scatter(cands_adaptive_NCC(1,2),cands_adaptive_NCC(1,1),80,"filled","square",'white');
% plot(cands_adaptive_NCC(:,2),cands_adaptive_NCC(:,1),'white+','MarkerSize',12,'LineWidth',1.5);
% scatter(H_first_point_index(2),H_first_point_index(1),180,"filled","pentagram",'r');
% title("TUKEY score map");
% 
% subplot(2,3,5)
% imagesc(tercomMap);  hold on; grid on; colormap('jet'); colorbar; set(gca,'YDir','normal');
% title("Map");
% 
% 
% muz = 3;
% subplot(3,3,4)
% imagesc(score_map_DTW);  xlim([1,2*idx_e+1]); ylim([1,2*idx_n+1]);hold on; grid on; colormap('jet'); colorbar; set(gca,'YDir','normal');
% scatter(cands_adaptive_DTW(1,2),cands_adaptive_DTW(1,1),80,"filled","square",'white');
% plot(cands_adaptive_DTW(:,2),cands_adaptive_DTW(:,1),'white+','MarkerSize',12,'LineWidth',1.5);
% 
% subplot(3,3,5)
% imagesc(score_map_MI);  xlim([1,2*idx_e+1]); ylim([1,2*idx_n+1]);hold on; grid on; colormap('jet'); colorbar; set(gca,'YDir','normal');
% scatter(cands_adaptive_MI(1,2),cands_adaptive_MI(1,1),80,"filled","square",'white');
% plot(cands_adaptive_MI(:,2),cands_adaptive_MI(:,1),'white+','MarkerSize',12,'LineWidth',1.5);
% 
% subplot(3,3,6)
% imagesc(score_map_ROC);  xlim([1,2*idx_e+1]); ylim([1,2*idx_n+1]);hold on; grid on; colormap('jet'); colorbar; set(gca,'YDir','normal');
% scatter(cands_adaptive_ROC(1,2),cands_adaptive_ROC(1,1),80,"filled","square",'white');
% plot(cands_adaptive_ROC(:,2),cands_adaptive_ROC(:,1),'white+','MarkerSize',12,'LineWidth',1.5);
% 
% subplot(3,3,7)
% imagesc(score_map_CMB);  xlim([1,2*idx_e+1]); ylim([1,2*idx_n+1]);hold on; grid on; colormap('jet'); colorbar; set(gca,'YDir','normal');
% scatter(cands_adaptive_CMB(1,2),cands_adaptive_CMB(1,1),80,"filled","square",'white');
% plot(cands_adaptive_CMB(:,2),cands_adaptive_CMB(:,1),'white+','MarkerSize',12,'LineWidth',1.5);
% 





end



