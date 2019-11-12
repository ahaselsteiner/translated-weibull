pdTrue = TranslatedWeibull(2.776, 1.471, 0.8888);

x = [0:0.01:18];
f = pdTrue.pdf(x);

fig1 = figure('position', [100 100 450 280])
plot(x, f);
message = sprintf(['In Vanem and Bitner-Gregersen (2012), Fig. 14 \n' ...
     '(doi: 10.1016/j.apor.2012.05.006) \n' ...
     'the PDF peaks at ~2 m with density of ~0.26.']);
text(2, 0.3, message, 'horizontalalignment', ...
    'left', 'verticalalignment', 'bottom', 'fontsize', 8);
ylabel('Density (-)');
xlabel('Significant wave height (m)');
ylim([0 0.4]);
box off


n = 1000;
nOfSamples = 20;

alphaEstimated = nan(nOfSamples, 1);
betaEstimated = nan(nOfSamples, 1);
gammaEstimated = nan(nOfSamples, 1);
for i = 1:nOfSamples
    sample = pdTrue.drawSample(n);
    pdEstimated(i) = TranslatedWeibull();
    pdEstimated(i).fitDist(sample);
    alphaEstimated(i) = pdEstimated(i).Alpha;
    betaEstimated(i) = pdEstimated(i).Beta;
    gammaEstimated(i) = pdEstimated(i).Gamma;
end

fig2 = figure('position', [100 100 400, 130]);
subplot(1, 3, 1)
hold on
plot([0.5 1.5], [1 1], '-k')
boxplot(alphaEstimated, {'$$\hat{\alpha}$$'})
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
text(1.15, 1, [num2str(mean(alphaEstimated), '%1.3f') '+-' ...
    num2str(std(alphaEstimated), '%1.3f')], 'fontsize', 8, ...
    'verticalalignment', 'bottom'); 
box off

subplot(1, 3, 2)
hold on
plot([0.5 1.5], [1 1], '-k')
boxplot(betaEstimated, {'$$\hat{\beta}$$'})
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
text(1.15, 1, [num2str(mean(betaEstimated), '%1.3f') '+-' ...
    num2str(std(betaEstimated), '%1.3f')], 'fontsize', 8, ...
    'verticalalignment', 'bottom');
box off

subplot(1, 3, 3)
hold on
plot([0.5 1.5], [2 2], '-k')
boxplot(gammaEstimated, {'$$\hat{\gamma}$$'})
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
text(1.15, 2, [num2str(mean(gammaEstimated), '%1.3f') '+-' ...
    num2str(std(gammaEstimated), '%1.3f')], 'fontsize', 8, ...
    'verticalalignment', 'bottom'); 
box off
