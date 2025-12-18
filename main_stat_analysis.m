close all;clear;clc; format longG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rt = load('feature_data.mat');

results_table = rmmissing(rt.results_table);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yMAD_all_N = results_table.absolute_error_MAD;

sc = unique(results_table.sampleCount);
meanErr = arrayfun( @(n) mean(yMAD_all_N(results_table.sampleCount == n)),sc);
medianErr = arrayfun( @(n) median(yMAD_all_N(results_table.sampleCount == n)),sc);
stdErr = arrayfun( @(n) std(yMAD_all_N(results_table.sampleCount == n)),sc);

figure()
errorbar(sc,meanErr,stdErr,'o-','LineWidth',2,'MarkerSize',10,'CapSize',12); grid on;hold on;
xlabel('Sample Count (N)',"FontSize",13);
ylabel('Error (m)',"FontSize",13);
legend("Mean Error",'FontSize',12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sample count set as 40
results_table = results_table(results_table.sampleCount ==40,:);
results_table.locationIdx= [];
results_table.sampleCount = [];


% Error metrics extracted from table,
% Now table includes only features
yMAD = results_table.absolute_error_MAD;
yMSD = results_table.absolute_error_MSD;
yNCC = results_table.absolute_error_NCC;
yDTW = results_table.absolute_error_DTW;
yMI  = results_table.absolute_error_MI;
yROC = results_table.absolute_error_ROC;
yPRO = results_table.absolute_error_PRO;

results_table.absolute_error_MAD = [];
results_table.absolute_error_MSD = [];
results_table.absolute_error_NCC = [];
results_table.absolute_error_DTW = [];
results_table.absolute_error_MI  = [];
results_table.absolute_error_ROC = [];
results_table.absolute_error_PRO = [];


yGenel = [yMAD, yMSD, yNCC, yDTW,  yMI , yROC];
yGenelTitles= ["MAD", "MSD", "NCC", "DTW", "MI" , "ROC"];


figure()
sgtitle("Error Distributions of Matching Metrics","FontSize",15)
for i= 1:6
    subplot(2,3,i)

    histogram(yGenel(:,i),"Normalization","probability" ,FaceColor=[0.2 0.4 0.7]); xlabel("Absolute Error (m)", "FontSize",13), ylabel("Occurance Density", "FontSize",13)
    title(sprintf("%s",yGenelTitles(i)),"FontSize",15);
    text(450.7,0.55,sprintf("mean:    %.2f\nmedian: %.2f",mean(yGenel(:,i)), median(yGenel(:,i)) ) ,"FontSize",12)
    xlim([0 max(yGenel(:,i))]);
    ylim([0 0.65])
end


figure()
subplot(1,3,1)
histogram(yMAD,"Normalization","probability" ,FaceColor=[0.2 0.4 0.7]); 
title("MAD","FontSize",17);
xlabel("Absolute Error (m)", "FontSize",15); ylabel("Occurance Density", "FontSize",15);
text(450.7,0.45,sprintf("mean:    %.2f\nmedian: %.2f",mean(yMAD), median(yMAD) ) ,"FontSize",12)
xlim([0 max(yMAD)]);    ylim([0 0.5])

subplot(1,3,2)
histogram(yDTW,"Normalization","probability" ,FaceColor=[0.2 0.4 0.7]); 
title("DTW","FontSize",17);
xlabel("Absolute Error (m)", "FontSize",15); ylabel("Occurance Density", "FontSize",15);
text(450.7,0.45,sprintf("mean:    %.2f\nmedian: %.2f",mean(yDTW), median(yDTW) ) ,"FontSize",12)
xlim([0 max(yDTW)]);    ylim([0 0.5])

subplot(1,3,3)
histogram(yPRO,"Normalization","probability" ,FaceColor=[0.2 0.4 0.7]); 
title("PRO","FontSize",17);
xlabel("Absolute Error (m)", "FontSize",15); ylabel("Occurance Density", "FontSize",15);
text(450.7,0.45,sprintf("mean:    %.2f\nmedian: %.2f",mean(yPRO), median(yPRO) ) ,"FontSize",12)
xlim([0 max(yPRO)]);    ylim([0 0.5])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ratioMAD = results_table.ratioMAD;
ratioDTW = results_table.ratioDTW;

% To simplfy only ratioMAD holded in the features.
results_table.ratioMSD = [];
results_table.ratioNCC = [];
results_table.ratioDTW = [];
results_table.ratioMI = [];
results_table.ratioROC = [];

results_table = rmmissing(results_table);
% X = results_table{:, vartype("numeric")};

X = table2array(results_table);
Xz = zscore(X);

featureNames = results_table.Properties.VariableNames;
lowerTh = 90;
upperTh = lowerTh;
numFeatures = size(featureNames,2);


%% X-X, feature korelasyonlar
R_p = corr(Xz, 'Type', 'Pearson', 'Rows', 'complete');
R_s = corr(Xz, 'Type', 'Spearman', 'Rows', 'complete');

XvarNames = results_table.Properties.VariableNames;
figure()
hm_p = heatmap(XvarNames,XvarNames,R_p, ...
    'ColorMap', parula, ...
    'ColorLimits',[-1 1], ...
    'MissingDataColor', [0.8 0.8 0.8]);
hm_p.ColorbarVisible = 'on';
hm_p.Title = 'Pearson Korelasyon Matrisi';

figure()
hm_s = heatmap(XvarNames,XvarNames,R_s, ...
    'ColorMap', parula, ...
    'ColorLimits',[-1 1], ...
    'MissingDataColor', [0.8 0.8 0.8]);
hm_s.ColorbarVisible = 'on';
hm_s.Title = 'Spearman Korelasyon Matrisi';


thLow = 0.80;
pairs = {};
elim = {};

for i = 1:size(R_p,1)
    for j = 1:size(R_s,1)
        if i == j
            continue
        end
        if abs(R_p(i,j)) >= thLow || abs(R_s(i,j)) >= thLow

            pairs(end+1,:) = {XvarNames{i},XvarNames{j}, R_p(i,j), R_s(i,j)};

            std_i = std(Xz(:,i),'omitnan');
            std_j = std(Xz(:,j),'omitnan');

            if std_i < std_j
                elim{end+1} = XvarNames{i};
            else
                elim{end+1} = XvarNames{j};
            end


        end
    end
end

corrTable = cell2table(pairs, ...
    "VariableNames",{'Feature 1', 'Feature 2','Rp','Rs'});
elim = unique(elim)';
disp('Korelasyondan Elenecek featurelar');
disp(elim);

results_table.avg_score_MAD = [];
results_table.mapEntropy = [];
results_table.mapFFTMaxPowFs = [];
results_table.mapFFTMeanPowFs= [];
results_table.mapFFTstdPowFs = [];
results_table.mapIqr = [];
results_table.mapMean = [];
results_table.mapMeanSlope = [];
results_table.mapMin = [];
results_table.mapRange = [];
results_table.mapTPImean = [];
results_table.mapTRImean = [];
results_table.mapVariance = [];
results_table.scoreMapFFTMaxPowFs = [];
results_table.scoreMapFFTMeanPowFs= [];
results_table.scoreMapFFTMedianPowFs= [];
results_table.scoreMapFFTstdPowFs = [];
results_table.scoreMapIqr = [];
results_table.scoreMapMax = [];
results_table.scoreMapMean = [];
results_table.scoreMapMeanPlanaryCurvature = [];
results_table.scoreMapMeanProfileCurvature = [];
results_table.scoreMapMeanSlope = [];
results_table.scoreMapMin = [];
results_table.scoreMapRange = [];
results_table.scoreMapSkewness = [];
results_table.scoreMapVariance = [];



X = table2array(results_table);
Xz = zscore(X);

featureNames = results_table.Properties.VariableNames;
fprintf("\n After pairvise correlation elemination %d features remained.\n", size(featureNames,2));

%% Top iin MI
p = numel(featureNames);

nbins = round(sqrt(size(Xz,1)));
MI_top = zeros(p,1);
for j = 1:p
    MI_top(j) = mutual_info(Xz(:,j), yMAD, nbins);
end

T_MI = table(featureNames', MI_top, 'VariableNames',{'Feature','MI'});
disp(T_MI);

MI_norm = (MI_top - min(MI_top)) ./ (max(MI_top) - min(MI_top));
MI_th = 0.25;
elimIdxMI = MI_norm <MI_th;

disp('MIa gre elenecek featurelar');
disp(featureNames(elimIdxMI)');

results_table.mapTRIARSclassIndex = [];
results_table.mapTPIPVS = [];
results_table.mapTPISmax = [];
results_table.mapSkewness= [];
results_table.mapKurtosis = [];
results_table.mapFFTMedianPowFs = [];
results_table.scoreMapKurtosis = [];


X = table2array(results_table);
Xz = zscore(X);

featureNames = results_table.Properties.VariableNames;

fprintf("\n After MI relevance elemination %d features remained.\n", size(featureNames,2));

%% AUC
yClass = double(yMAD < 90);

p = size(Xz,2);
AUCvals = zeros(p,1);

for j = 1:p
    [~,~,~, AUC] = perfcurve(yClass,Xz(:,j),1);
    AUCvals(j) = AUC;
end

T_AUC = table(featureNames(:), AUCvals, ...
    'VariableNames', {'Feature', 'AUC'});

T_AUC = sortrows(T_AUC,'AUC','descend');


T_AUC_elimIdx = T_AUC.AUC < 0.5;
T_AUC_elimFeatures = T_AUC.Feature(T_AUC_elimIdx);

disp('Feature AUC sralamas:');
disp(T_AUC);

disp('AUC dan elenecek Featurelar');
disp(T_AUC_elimFeatures);

results_table.mapMeanProfileCurvature = [];
results_table.mapMeanPlanaryCurvature = [];
results_table.mapFFTRadialSpectralEntropy = [];
results_table.scoreMapFFTSpectralEntrpoy = [];
results_table.scoreMapFFTbw_u = [];

X = table2array(results_table);
Xz = zscore(X);

featureNames = results_table.Properties.VariableNames;

fprintf("\n After AUC relevance elemination %d features remained.\n", size(featureNames,2));


%% mRMR
[mRMR_idx,mRMR_scores] = fscmrmr(X,yClass);
featureNames(mRMR_idx)
mRMR_scores(mRMR_idx);

mRMR_elimIdx = mRMR_scores < 0.02;
mRMR_elimFeatures = T_AUC.Feature(mRMR_elimIdx)


results_table.mapMax = [];
results_table.mapTPISAR = [];

X = table2array(results_table);
Xz = zscore(X);
featureNames = results_table.Properties.VariableNames;

fprintf("\n After nRMR elemination %d features remained.\n", size(featureNames,2));


%% End of Feature Screening
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
results_table.ratioMAD = ratioMAD;
results_table.ratioDTW = ratioDTW;

ErrorTh = 90;

y_MAD_idx = yMAD <= ErrorTh;
y_DTW_idx = yDTW <= ErrorTh;

X_list = {ratioMAD, ratioDTW, results_table.roughness, results_table.roughness, results_table.roughnessDiff, results_table.roughnessDiff };
Y_list = {y_MAD_idx, y_DTW_idx, y_MAD_idx, y_DTW_idx, y_MAD_idx, y_DTW_idx};
names = {'ratioMAD - MAD absolute error', ...
         'ratioDTW - DTW absolute error', ...
         'roughness - MAD absolute error', ...
         'roughness- DTW absolute error', ...
         'roughness difference - MAD absolute error', ...
         'roughness difference - DTW absolute error'};

th_best = zeros(6,1);
J_best = zeros(6,1);
TPR_best = zeros(6,1);
FPR_best = zeros(6,1);

figure();
for k = 1:6

    xvals = X_list{k};
    y     = Y_list{k};
    
    tmin = min(xvals);
    tmax = max(xvals);
    thresholds = linspace(tmin,tmax,500);
    
    TPR_rr = zeros(size(thresholds));
    FPR_rr = zeros(size(thresholds));
    J_rr   = zeros(size(thresholds));
    Precision_rr = zeros(size(thresholds));
    Recall_rr = zeros(size(thresholds));
    Accuracy_rr = zeros(size(thresholds));

    for i = 1:length(thresholds)
    
        th = thresholds(i);

        preds = xvals > th;
%         N = numel(xvals);
        TP = sum(preds == 1 & y == 1);
        FP = sum(preds == 1 & y == 0);
        FN = sum(preds == 0 & y == 1);
        TN = sum(preds == 0 & y == 0);

%         Prec(i) = (TP+FP) / (TP + FP + FN + TN + eps);
%         Rec(i) = TP / (TP + FN + eps);
        Precision_rr(i) = (TP) / (TP + FP + eps);
        Recall_rr(i) = (TP) / (TP + FN + eps);
%         Accuracy(i) = ( TP + TN) / (TP + FP + FN + TN + eps);
        Accuracy_rr(i) = (FP) / (TP + FP + eps);

        rateTP(i) = (TP) / (TP + FP + FN + TN + eps);
        rateFP(i) = (FP) / (TP + FP + FN + TN + eps);
        rateFN(i) = (FN) / (TP + FP + FN + TN + eps);
        rateTN(i) = (TN) / (TP + FP + FN + TN + eps);

        TPR_rr(i) = TP / ( TP + FN + eps);
        FPR_rr(i) = FP / ( FP + TN + eps);    
        J_rr(i) = TPR_rr(i) - FPR_rr(i);
%         Keep(i) = (TP+TN) /N;
     end

    REC_TARGET = 0.90;
    mRMR_idx = find(Recall_rr <= REC_TARGET,1,'first');


    [Jmax, ggggg] = max(J_rr);
    th_star = thresholds(mRMR_idx);

    th_best(k) = th_star;
    J_best(k) = Jmax;
    TPR_best(k) = TPR_rr(mRMR_idx);
    FPR_best(k) = FPR_rr(mRMR_idx);

    subplot(3,2,k)
%     plot(thresholds,rateTP, '-', 'LineWidth',2); hold on;
%     plot(thresholds,rateFP, '-', 'LineWidth',2); hold on;
% %     plot(thresholds,rateFN, '-', 'LineWidth',2); hold on;
% %     plot(thresholds,rateTN, '-', 'LineWidth',2); hold on;
%     legend("rateTP","rateFP","rateFN","rateTN");


    plot(thresholds,Precision_rr, '-', 'LineWidth',2); hold on;
    plot(thresholds,Recall_rr, '-', 'LineWidth',2); hold on;
    plot(thresholds,Accuracy_rr, '-', 'LineWidth',2); hold on;
%     plot(FPR(idx),TPR(idx), 'ro', 'MarkerSize',8,'LineWidth',2);
    legend("Precision","Recall","Accuracy");
    grid on; %axis square;
    xlabel('feature degeri'); ylabel('Oran');
    title(sprintf('%s\nJ=%.2f, TH=%.3g', names{k}, Jmax, th_star),'Interpreter','none');
end

%%
% TH_ratioMAD = 18.9;
% TH_ratioDTW = 17;
% TH_roughnessMAD = 21.8;
% TH_roughnessDTW = 23.9;
% TH_roughnessDiffMAD = 3.92;
% TH_roughnessDiffDTW = 4.33;

TH_ratioMAD = th_best(1);
TH_ratioDTW = th_best(2);
TH_roughnessMAD = th_best(3);
TH_roughnessDTW = th_best(5);
TH_roughnessDiffMAD = th_best(5);
TH_roughnessDiffDTW = th_best(6);

abstain = zeros(size(results_table.ratioMAD,1),1);
HATA = zeros(size(results_table.ratioMAD,1),1);
for i = 1:size(results_table.ratioMAD,1)

    chooseMAD = 0;
    chooseDTW = 0;
    if results_table.ratioMAD > results_table.ratioDTW 
        ratioPROO = results_table.ratioMAD;
        chooseMAD = 1;
    else
        ratioPROO = results_table.ratioDTW;
        chooseDTW = 1;
    end

    if chooseMAD == 1
        a = results_table.ratioMAD(i) < TH_ratioMAD;
        b = results_table.roughness(i) < TH_roughnessMAD;
        c = results_table.roughnessDiff(i) < TH_roughnessDiffMAD;
        abstain(i) = a | b | c;
        if abstain(i)
            HATA(i) = 0;
        else
            HATA(i) = yMAD(i);
        end
    else
        a = results_table.ratioDTW(i) < TH_ratioDTW;
        b = results_table.roughness(i) < TH_roughnessDTW;
        c = results_table.roughnessDiff(i) < TH_roughnessDiffDTW;
        abstain(i) = a | b | c;
        if abstain(i)
            HATA(i) = 0;
        else
            HATA(i) = yDTW(i);
        end
    end

end

figure()
histogram(HATA); hold on;
histogram(yMAD); hold on;
histogram(yDTW); hold on;

sum(HATA > 90)
sum(yMAD > 90)
sum(yDTW > 90)

figure()
histogram(abstain,2,"Normalization","probability")



%% Roughness && Roughnes Difference plots

X = table2array(results_table);
Xz = zscore(X);

featureNames = results_table.Properties.VariableNames;
numFeatures = size(featureNames,2);


X_list = {results_table.ratioMAD, results_table.ratioDTW, results_table.roughness, results_table.roughness, results_table.roughnessDiff, results_table.roughnessDiff, results_table.scoreMapEntropy,results_table.scoreMapEntropy };
% Y_list = {y_MAD_idx, y_DTW_idx, y_MAD_idx, y_DTW_idx, y_MAD_idx, y_DTW_idx};
Y_list = {yMAD, yDTW, yMAD, yDTW, yMAD, yDTW, yMAD, yDTW};
names = {'ratioMAD - MAD absolute error', ...
         'ratioDTW - DTW absolute error', ...
         'roughness - MAD absolute error', ...
         'roughness- DTW absolute error', ...
         'roughness difference - MAD absolute error', ...
         'roughness difference - DTW absolute error', ...
         'scEntropy - MAD absolute error', ...
         'scEntropy - DTW absolute error'};


figure();
for k = 1:8

    xvals = X_list{k};
    y     = Y_list{k};

    edges = linspace(min(xvals(:)),max(xvals(:)),151);
    
    lowErr = y <90;
    highErr = ~lowErr;
    [countLow,~] = histcounts(xvals(lowErr),edges);
    [countHigh,~] = histcounts(xvals(highErr ),edges);
    N = size(xvals,1);
    propLow = countLow/N;
    propHigh = countHigh/N;
    subplot(4,2,k)
    bar(edges(1:end-1), [propLow' propHigh'], 'stacked');
    title(sprintf('Feature %d : %s',k,names{k}),"FontSize",15)
    ylabel(" oran, tum veri");
    legend('< 90', '>= 90');


end



%% RatioMAD & RadioDTW plots

X = table2array(results_table);
Xz = zscore(X);

featureNames = results_table.Properties.VariableNames;
numFeatures = size(featureNames,2);


X_list = {results_table.ratioMAD, results_table.ratioDTW, results_table.roughness, results_table.roughness, results_table.roughnessDiff, results_table.roughnessDiff };
% Y_list = {y_MAD_idx, y_DTW_idx, y_MAD_idx, y_DTW_idx, y_MAD_idx, y_DTW_idx};
Y_list = {yMAD, yDTW, yMAD, yDTW, yMAD, yDTW};
names = {'ratioMAD - MAD absolute error', ...
         'ratioDTW - DTW absolute error', ...
         'roughness - MAD absolute error', ...
         'roughness- DTW absolute error', ...
         'roughness difference - MAD absolute error', ...
         'roughness difference - DTW absolute error'};


figure();
for k = 1:2

    xvals = X_list{k};
    y     = Y_list{k};

    edges = linspace(min(xvals(:)),max(xvals(:)),101);
    
    lowErr = y <90;
    highErr = ~lowErr;
    [countLow,~] = histcounts(xvals(lowErr),edges);
    [countHigh,~] = histcounts(xvals(highErr ),edges);
    N = size(xvals,1);
    propLow = countLow/N;
    propHigh = countHigh/N;

    subplot(1,2,k)
    barHead = bar(edges(1:end-1), [propLow' propHigh'], 'stacked');
    barHead(2).FaceColor =[0.8 0 0];

    title("Minima Ratio - Error Histogram","FontSize",15);% sprintf('%s',names{k}),"FontSize",15)
    ylabel("Occurance Density");  
    if k == 1 
        xlabel("Minima Ratio (MAD)"); 
        title("Minima Ratio - Error Histogram (MAD)","FontSize",12);
    else 
        xlabel("Minima Ratio (DTW)");
        title("Minima Ratio - Error Histogram (DTW)","FontSize",12);
    end
    
    cumLow = cumsum(countLow);
    cumHigh = cumsum(countHigh);
    ratioCumSum = cumHigh ./ cumLow ;
    
    yyaxis right
    plot(edges(1:end-1),ratioCumSum, "LineWidth",2,"Color",[0, 0.8, 0]);
    ax = gca;
    ax.YColor = [0, 0.5, 0];
    ylabel("Cumulative Ratio");
    legend('Absolute Error < 90 m', 'Absolute Error >= 90 m ','Cumulative Ratio of Bad over Good Solutions','FontSize',9);

end


%%
X = table2array(results_table);
Xz = zscore(X);

featureNames = results_table.Properties.VariableNames;
numFeatures = size(featureNames,2);


X_list = {results_table.scoreMapEntropy, results_table.scoreMapEntropy, results_table.roughness, results_table.roughness, results_table.roughnessDiff, results_table.roughnessDiff };
% Y_list = {y_MAD_idx, y_DTW_idx, y_MAD_idx, y_DTW_idx, y_MAD_idx, y_DTW_idx};
Y_list = {yMAD, yDTW, yMAD, yDTW, yMAD, yDTW};
names = {'ratioMAD - MAD absolute error', ...
         'ratioDTW - DTW absolute error', ...
         'roughness - MAD absolute error', ...
         'roughness- DTW absolute error', ...
         'roughness difference - MAD absolute error', ...
         'roughness difference - DTW absolute error'};


figure();
for k = 1:2

    xvals = X_list{k};
    y     = Y_list{k};

    edges = linspace(min(xvals(:)),max(xvals(:)),101);
    
    lowErr = y <90;
    highErr = ~lowErr;
    [countLow,~] = histcounts(xvals(lowErr),edges);
    [countHigh,~] = histcounts(xvals(highErr ),edges);
    N = size(xvals,1);
    propLow = countLow/N;
    propHigh = countHigh/N;

    subplot(2,2,k*2-1)
    barHead = bar(edges(1:end-1), [propLow' propHigh'], 'stacked');
    barHead(2).FaceColor =[0.8 0 0];

    title("Minima Ratio - Error Histogram","FontSize",15);% sprintf('%s',names{k}),"FontSize",15)
    ylabel("Occurance Density");  
    if k == 1 
        xlabel("Minima Ratio (MAD)"); 
        title("scent - Error Histogram (MAD)","FontSize",12);
    else 
        xlabel("Minima Ratio (DTW)");
        title("scent - Error Histogram (DTW)","FontSize",12);
    end
    
    cumLow = cumsum(countLow);
    cumHigh = cumsum(countHigh);
    ratioCumSum = cumHigh ./ cumLow ;
    
    yyaxis right
    plot(edges(1:end-1),ratioCumSum, "LineWidth",2,"Color",[0, 0.8, 0]);
    ax = gca;
    ax.YColor = [0, 0.5, 0];
    ylabel("Cumulative Ratio");
    legend('Absolute Error < 90 m', 'Absolute Error >= 90 m ','Cumulative Ratio of Bad over Good Solutions','FontSize',9);

end

sum(results_table.scoreMapEntropy < 6)

for k = 1:2

    xvals = X_list{k};
    y     = Y_list{k};

    edges = linspace(4,max(xvals(:)),101);
    
    lowErr = y <90;
    highErr = ~lowErr;
    [countLow,~] = histcounts(xvals(lowErr),edges);
    [countHigh,~] = histcounts(xvals(highErr ),edges);
    N = size(xvals,1);
    propLow = countLow/N;
    propHigh = countHigh/N;

    subplot(2,2,2*k)
    barHead = bar(edges(1:end-1), [propLow' propHigh'], 'stacked');
    barHead(2).FaceColor =[0.8 0 0];

    title("Minima Ratio - Error Histogram","FontSize",15);% sprintf('%s',names{k}),"FontSize",15)
    ylabel("Occurance Density");  
    if k == 1 
        xlabel("Minima Ratio (MAD)"); 
        title("scent - Error Histogram (MAD)","FontSize",12);
    else 
        xlabel("Minima Ratio (DTW)");
        title("scent - Error Histogram (DTW)","FontSize",12);
    end
    
    cumLow = cumsum(countLow);
    cumHigh = cumsum(countHigh);
    ratioCumSum = cumHigh ./ cumLow ;
    
    yyaxis right
    plot(edges(1:end-1),ratioCumSum, "LineWidth",2,"Color",[0, 0.8, 0]);
    ax = gca;
    ax.YColor = [0, 0.5, 0];
    ylabel("Cumulative Ratio");
    legend('Absolute Error < 90 m', 'Absolute Error >= 90 m ','Cumulative Ratio of Bad over Good Solutions','FontSize',9);

end


%%
X = table2array(results_table);
Xz = zscore(X);

featureNames = results_table.Properties.VariableNames;
numFeatures = size(featureNames,2);


X_list = {results_table.roughness, results_table.roughness, results_table.roughnessDiff, results_table.roughnessDiff };
% Y_list = {y_MAD_idx, y_DTW_idx, y_MAD_idx, y_DTW_idx, y_MAD_idx, y_DTW_idx};
Y_list = {yMAD, yDTW, yMAD, yDTW};
names = {'roughness - MAD absolute error', ...
         'roughness- DTW absolute error', ...
         'roughness difference - MAD absolute error', ...
         'roughness difference - DTW absolute error'};

TableNames = {"Roughness - Error Histogram (MAD)", ...
              "Roughness - Error Histogram (DTW)", ...
              "Roughness Difference - Error Histogram (MAD)", ...
              "Roughness Difference - Error Histogram (DTW)"};

figure();
for k = 1:4

    xvals = X_list{k};
    y     = Y_list{k};

    edges = linspace(min(xvals(:)),max(xvals(:)),101);
    
    lowErr = y <90;
    highErr = ~lowErr;
    [countLow,~] = histcounts(xvals(lowErr),edges);
    [countHigh,~] = histcounts(xvals(highErr ),edges);
    N = size(xvals,1);
    propLow = countLow/N;
    propHigh = countHigh/N;

    subplot(4,2,k*2-1)
    barHead = bar(edges(1:end-1), [propLow' propHigh'], 'stacked');
    barHead(2).FaceColor =[0.8 0 0];

    title(sprintf('%s',TableNames{k}),"FontSize",12)
    ylabel("Occurance Density");  
    % if k == 1 
    %     xlabel("Minima Ratio (MAD)"); 
    %     title("Minima Ratio - Error Histogram (MAD)","FontSize",12);
    % else 
    %     xlabel("Minima Ratio (DTW)");
    %     title("Minima Ratio - Error Histogram (DTW)","FontSize",12);
    % end
    
    cumLow = cumsum(countLow);
    cumHigh = cumsum(countHigh);
    ratioCumSum = cumHigh ./ cumLow ;
    
    yyaxis right
    plot(edges(1:end-1),ratioCumSum, "LineWidth",2,"Color",[0, 0.8, 0]);
    ax = gca;
    ax.YColor = [0, 0.5, 0];
    ylabel("Cumulative Ratio");
    legend('Absolute Error < 90 m', 'Absolute Error >= 90 m ','Cumulative Ratio of Bad over Good Solutions','FontSize',8);

end


for k = 1:4

    xvals = X_list{k};
    y     = Y_list{k};

    edges = linspace(min(xvals(:)),max(xvals(:))/400,101);
    
    lowErr = y <90;
    highErr = ~lowErr;
    [countLow,~] = histcounts(xvals(lowErr),edges);
    [countHigh,~] = histcounts(xvals(highErr ),edges);
    N = size(xvals,1);
    propLow = countLow/N;
    propHigh = countHigh/N;

    subplot(4,2,k*2)
    barHead = bar(edges(1:end-1), [propLow' propHigh'], 'stacked');
    barHead(2).FaceColor =[0.8 0 0];

    title(sprintf('Zoomed %s',TableNames{k}),"FontSize",12)
    ylabel("Occurance Density");  
    % if k == 1 
    %     xlabel("Minima Ratio (MAD)"); 
    %     title("Minima Ratio - Error Histogram (MAD)","FontSize",12);
    % else 
    %     xlabel("Minima Ratio (DTW)");
    %     title("Minima Ratio - Error Histogram (DTW)","FontSize",12);
    % end
    
    cumLow = cumsum(countLow);
    cumHigh = cumsum(countHigh);
    ratioCumSum = cumHigh ./ cumLow ;
    

    
    yyaxis right
    plot(edges(1:end-1),ratioCumSum, "LineWidth",2,"Color",[0, 0.8, 0]); hold on; %grid on;
    ax = gca;
    ax.YColor = [0, 0.5, 0];
        targets = [5 2 1];
    for t = targets
        [~, idx_target] = min(abs(ratioCumSum - t));
        plot(edges(idx_target), ratioCumSum(idx_target),'white*');
        xline(edges(idx_target),'black--','LineWidth',1.2);
    end
    ylabel("Cumulative Ratio");
    legend('Absolute Error < 90 m', 'Absolute Error >= 90 m ','Cumulative Ratio of Bad over Good Solutions','FontSize',8);

end



%%

TH_ratioMAD         = linspace(0.1,100,101);
TH_roughnessMAD     = linspace(0,5,101);
TH_roughnessDiffMAD = linspace(0,1,101);
TH_ratioDTW         = linspace(0.1,100,101);
TH_roughnessDTW     = linspace(0,5,101);
TH_roughnessDiffDTW = linspace(0,1,101);

for i = 1:101
    for j = 1:size(results_table.roughness,1)

        chooseMAD = 0;
        chooseDTW = 0;
        if results_table.ratioMAD(j)> results_table.ratioDTW(j)
            a = results_table.ratioMAD(j) < TH_ratioMAD(i);
            b = results_table.roughness(j) < TH_roughnessMAD(i);
            c = results_table.roughnessDiff(j) < TH_roughnessDiffMAD(i);

            abstain(j) = a | b | c;
            errrrr(j) = yMAD(j);
            
        else
            a = results_table.ratioDTW(j) < TH_ratioDTW(i);
            b = results_table.roughness(j) < TH_roughnessDTW(i);
            c = results_table.roughnessDiff(j) < TH_roughnessDiffDTW(i);
            abstain(j) = a | b | c;
            errrrr(j) = yDTW(j);
            
        end
        
        
        
    end
    go = 0;

    flag_dusukHata = errrrr < 90;
    [counts,~,~,labelsResult] = crosstab(flag_dusukHata', abstain');
    countsNorm = counts./sum(counts(:));
    TP = countsNorm(2,1);
    FP = countsNorm(1,1);
    FN = countsNorm(2,2);
    TN = countsNorm(1,2);
    precision(i)   = TP / (TP + FP);
    recall(i)      = TP / (TP + FN); 
    accuracy(i)    = (TP + TN) /(TP + TN + FP +FN);
    falseFix(i)    = FP;

    TP_count = counts(2,1);
    FP_count = counts(1,1);
    FN_count = counts(2,2);
    TN_count = counts(1,2);

    TPR(i) = TP_count / (TP_count + FN_count + eps);
    FPR(i) = FP_count / (FP_count + TN_count + eps);
end

figure()
plot(FPR, TPR,'-o','Linewidth',1.5);hold on;
plot([0 1], [0 1], '--k');
grid on;
xlabel('False Positive Rate','FontSize',12); ylabel('True Positive Rate','FontSize',12);
title('ROC Curve for B/G Threshold Sweep','FontSize',15);

[~, idxFPR] = sort(FPR);
AUC = trapz(FPR(idxFPR), TPR(idxFPR));
fprintf('ROC AUC = %.3f\n', AUC);

highErrhlight_idx = [6 18 34];
BG= [5,2,1];
counter = 0;
for id = highErrhlight_idx
    counter = counter +1;
    plot(FPR(id), TPR(id), 'rs', 'MarkerSize',12,'LineWidth',2);
    text(FPR(id), TPR(id)-0.015, sprintf(' B/G =%.1f',BG(counter)));
end

tpr_fpr_rate = diff(TPR)./diff(FPR);
figure()
plot(tpr_fpr_rate)

figure()
% subplot(1,2,1);
plot(precision,"LineWidth",2); hold on; grid on;
plot(accuracy,"LineWidth",2);
plot(recall,"LineWidth",2);
plot(falseFix,"LineWidth",2);
xline(6,'black--','LineWidth',2);
xline(18,'black:','LineWidth',2);
xline(34,'black-.','LineWidth',2);
legend("Precision","Accuracy","Recall","False Fix","B/G = 5","B/G = 2","B/G = 1","location","best");
xlabel("Iteration","FontSize",13); ylabel("Value","FontSize",13); xlim([1 100]);
title("Performance Metrics vs. Iteration","FontSize",15);


figure()
plot(recall,precision,"LineWidth",2); hold on; grid on;
plot(recall,accuracy,"LineWidth",2); hold on; grid on;
plot(recall,falseFix,"LineWidth",2); hold on; grid on;
xlabel("Recall","FontSize",13); ylabel("Value","FontSize",13);
xline(0.9945,'black--','LineWidth',2);
xline(0.9756,'black:','LineWidth',2);
xline(0.9075,'black-.','LineWidth',2);
legend("Precision - Recall","Accuracy - Recall","False Fix - Recall" ,"B/G = 5","B/G = 2","B/G = 1","location","best");
title("Performance Metrics vs. Recall","FontSize",15);




%% Bar plots of Confusion Matrix Elements

skipRatio = results_table.ratioDTW <= 5 | results_table.ratioMAD <= 4 |results_table.ratioMAD >= 98 |results_table.ratioDTW >= 99 ; 
skipRoughness =  results_table.roughness <= 0.5; % 1.3;            5 
skipRoughnessDiff =  results_table.roughnessDiff <= 0.2 ;

skipSolIdx = skipRatio | skipRoughness | skipRoughnessDiff; %| skipAvgDiffMAD;% | skipBestScoreMAD;


flag_dusukHata = yDTW < 90;
[counts,~,~,labelsResult] = crosstab(flag_dusukHata, skipSolIdx);

max(yDTW(~skipSolIdx & ~flag_dusukHata))
max(yDTW)

figure()
histogram(yDTW(~skipSolIdx & ~flag_dusukHata))

countsNorm = counts./sum(counts(:));

edges= [-0.5 0.5 1.5];
counts = histcounts2(flag_dusukHata,skipSolIdx,edges,edges);
FP_oran = counts(1,1)/sum(sum(counts));
TP_oran = counts(2,1)/sum(sum(counts));
[TP_oran, FP_oran];

figure()
subplot(1,3,1)
flagBar = bar(countsNorm,'grouped');
set(gca,"XTickLabel",["High Error (E>90)", "Low Error (E<90)"],'FontSize',8);
xlabel("Groundtruth Flag")
ylabel("Density")
ylim([0 1]);
title("B/G = 1")

legend("Fix", "Abstain",Location="northwest")

text(0.8,0.09, "FP","Color","black")
text(1.05,0.12,"TN","Color","black")
text(1.8,0.83, "TP","Color","black")
text(2.05,0.06,"FN","Color","black")



% 2. gorsel
skipRatio = results_table.ratioDTW <= 20 | results_table.ratioMAD <= 16 |results_table.ratioMAD >= 98 |results_table.ratioDTW >= 99 ; 
skipRoughness =  results_table.roughness <= 3; % 1.3;            5 
skipRoughnessDiff =  results_table.roughnessDiff <= 0.4 ;

skipSolIdx = skipRatio | skipRoughness | skipRoughnessDiff; %| skipAvgDiffMAD;% | skipBestScoreMAD;


flag_dusukHata = yDTW < 90;
[counts,~,~,labelsResult] = crosstab(flag_dusukHata, skipSolIdx);

max(yDTW(~skipSolIdx & ~flag_dusukHata))
max(yDTW)



countsNorm = counts./sum(counts(:));


edges= [-0.5 0.5 1.5];
counts = histcounts2(flag_dusukHata,skipSolIdx,edges,edges);
FP_oran = counts(1,1)/sum(sum(counts));
TP_oran = counts(2,1)/sum(sum(counts));


% figure()
subplot(1,3,2)
flagBar = bar(countsNorm,'grouped');
set(gca,"XTickLabel",["High Error (E>90)", "Low Error (E<90)"],'FontSize',8);
xlabel("Groundtruth Flag")
ylabel("Density")
legend("Fix", "Abstain",Location="northwest")
ylim([0 1]);
text(0.77,0.065, "FP","Color","black")
text(1.05,0.155,"TN","Color","black")
text(1.8,0.76, "TP","Color","black")
text(2.05,0.14,"FN","Color","black")
title("B/G = 1")


% 3. gorsel
skipRatio = results_table.ratioDTW <= 33 | results_table.ratioMAD <= 31 |results_table.ratioMAD >= 98 |results_table.ratioDTW >= 99 ; 
skipRoughness =  results_table.roughness <= 3.5; % 1.3;            5 
skipRoughnessDiff =  results_table.roughnessDiff <= 0.7 ;

skipSolIdx = skipRatio | skipRoughness | skipRoughnessDiff; %| skipAvgDiffMAD;% | skipBestScoreMAD;


flag_dusukHata = yDTW < 90;
[counts,~,~,labelsResult] = crosstab(flag_dusukHata, skipSolIdx);

max(yDTW(~skipSolIdx & ~flag_dusukHata));
max(yDTW);

countsNorm = counts./sum(counts(:));


edges= [-0.5 0.5 1.5];
counts = histcounts2(flag_dusukHata,skipSolIdx,edges,edges);
FP_oran = counts(1,1)/sum(sum(counts));
TP_oran = counts(2,1)/sum(sum(counts));

subplot(1,3,3)
flagBar = bar(countsNorm,'grouped');
set(gca,"XTickLabel",["High Error (E>90)", "Low Error (E<90)"],'FontSize',8);
xlabel("Groundtruth Flag")
ylabel("Density")
legend("Fix", "Abstain",Location="northwest")
ylim([0 1]);


text(0.8,0.05, "FP","Color","black")
text(1.05,0.16,"TN","Color","black")
text(1.8,0.65, "TP","Color","black")
text(2.05,0.23,"FN","Color","black")

title("B/G = 1")
sgtitle("Bar Plots of Confusion Matrix Elements")

%% 30m, 60m, 90m error, pozitive negative plots over ratio

for i = 1: numFeatures
    edges = linspace(min(X(:,i)),max(X(:,i)),51);

    figure()
    lowErr = yMAD <30;
    highErr = ~lowErr;
    [countLow,~] = histcounts(X(lowErr,i),edges);
    [countHigh,~] = histcounts(X(highErr ,i),edges);
    N = size(X,1);
    propLow = countLow/N;
    propHigh = countHigh/N;
    subplot(3,1,1)
    bar(edges(1:end-1), [propLow' propHigh'], 'stacked');
    xlabel(sprintf('Feature %d : %s',i,featureNames{i}),"FontSize",15)
    ylabel(" oran, tum veri");
    legend('< 30', '>= 30');

    lowErr = yMAD <60;
    highErr = ~lowErr;
    [countLow,~] = histcounts(X(lowErr,i),edges);
    [countHigh,~] = histcounts(X(highErr ,i),edges);
    N = size(X,1);
    propLow = countLow/N;
    propHigh = countHigh/N;
    subplot(3,1,2)
    bar(edges(1:end-1), [propLow' propHigh'], 'stacked');
    xlabel(sprintf('Feature %d : %s',i,featureNames{i}),"FontSize",15)
    ylabel(" oran, tum veri");
    legend('< 60', '>= 60');

    lowErr = yMAD <90;
    highErr = ~lowErr;
    [countLow,~] = histcounts(X(lowErr,i),edges);
    [countHigh,~] = histcounts(X(highErr ,i),edges);
    N = size(X,1);
    propLow = countLow/N;
    propHigh = countHigh/N;
    subplot(3,1,3)
    bar(edges(1:end-1), [propLow' propHigh'], 'stacked');
    xlabel(sprintf('Feature %d : %s',i,featureNames{i}),"FontSize",15)
    ylabel(" oran, tum veri");
    legend('< 90', '>= 90');
end


