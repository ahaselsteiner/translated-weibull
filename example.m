rng('default'); % For reproduceability.

pdTrue = TranslatedWeibull(2.8, 1.5, 0.9);
n = 10000;
sample = pdTrue.drawSample(n);

pdEstimated = TranslatedWeibull();
pdEstimated.fitDist(sample)
