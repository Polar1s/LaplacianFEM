% Data visualization for test.m
% by Beichen Li

% Global parameters
resPath = '../results';     % Path to mesh files
meshName = 'moomoo';    % Name of mesh
titleName = 'moomoo'; % Title of each figure
fileName = join([resPath,'/',meshName,'_res.mat']);

savePath = '../images';     % Path to images
legendFlag = 1;             % Switch for legend

% Load data (approximation error, time consumption)
load(fileName,'err','tc');

% Approximation error: display
figure('Position',[1000 600 350 300]);
plot(err([1 2 3 5 7 10],2:100)','LineWidth',1.5);
xlabel('Eigenfunction No.');
ylabel('Error');
if legendFlag == 1
    lgd = legend('n=1','n=2','n=3','n=5','n=7','n=10');
    lgd.Location = 'southeast';
end
set(gca,'YScale','log');
title(titleName);
% saveas(gcf,join([savePath,'/',meshName,'_err1.png']));

% Approximation error: comparison
ind = [8 32];
figure('Position',[1400 600 350 300]);
plot(reshape(err(:,ind),[10,4]),'LineWidth',2);
xlabel('n');
ylabel('Error');
if legendFlag == 1
    lgd = legend('No.8 H-O','No.8 REF','No.32 H-O','No.32 REF');
    lgd.Location = 'southwest';
end
title(titleName);
set(gca,'YScale','log');
% saveas(gcf,join([savePath,'/',meshName,'_err2.png']));

% Speed: comparison
figure('Position',[1000 200 350 200]);
bar(reshape(tc,[10 2]),'LineStyle','none');
xlabel('n');
ylabel('Time (s)');
if legendFlag == 1
    lgd = legend('H-O','REF');
    lgd.Location = 'northwest';
end
title(titleName);
% saveas(gcf,join([savePath,'/',meshName,'_speed.png']));

% Accuracy vs. speed
ind = 8;
figure('Position',[1400 200 250 200]);
plot(tc(1:10),err(1:10,ind),tc(11:20),err(11:20,ind),'LineWidth',2);
xlabel('Time (s)');
ylabel('Error');
if legendFlag == 1
    lgd = legend('H-O','REF');
end
title(titleName);
set(gca,'XScale','log');
set(gca,'YScale','log');
% saveas(gcf,join([savePath,'/',meshName,'_err_speed.png']));
