% This software was written for the publication
% "Predicting wave heights for marine design by prioritizing extreme 
% events in a global model" by A.F. Haselsteiner and K-D. Thoben, see
% https://arxiv.org/pdf/1911.12835.pdf .

classdef TranslatedWeibull < handle
% The translated Weibull distribution, sometimes also simply called 
% "3-parameter Weibull distribution" is a parametric probability
% distribution.
%
% We use the parameterization and variables names that are also used in 
% doi.org/10.1016/j.coastaleng.2017.03.002 .
   properties
      Alpha % Scale parameter.
      Beta % Shape parameter.
      Gamma % Location parameter.
      BootstrapParm % Parameters estimated using bootstrap.
      ParameterSE % Parameters' standard error estimated using bootstrapping.
      ParameterCI % Parameters'confidence interval standard error estimated using bootstrapping.
   end
   
   methods
      function obj = TranslatedWeibull(alpha, beta, gamma)
         if nargin > 2
            obj.Alpha = alpha;
            obj.Beta = beta;
            obj.Gamma = gamma;
         end
      end
      
      function parmHat = fitDist(this, sample)
          % Estimates the parameters of the distribution using maximum 
          % likelihood estimation.
          epsilon = 0.99999;
          max_location_param = epsilon * min(sample);
          start = [2 2 0];
          lb = [-Inf -Inf -Inf];
          ub = [Inf Inf max_location_param];
          parmHat = mle(sample, 'pdf', @(x, alpha, beta, gamma) ...
              this.pdf(sample, alpha, beta, gamma), ...
              'start', start, 'lower', lb, 'upper', ub);
          this.Alpha = parmHat(1);
          this.Beta = parmHat(2);
          this.Gamma = parmHat(3);
      end
      
      function [parmHat, pStd, pCi] = fitDistAndBootstrap(this, sample, B, alpha)
          % Estimates the parameters of the distribution using maximum 
          % likelihood estimation and estimate the parameters' uncertainty
          % using bootstrapping.
          %
          % For bootstrapping, see, for example, "An introduction to the
          % bootstrap" by B. Efron and R. J. Tibshirani (1993).
          if nargin < 4
              alpha = 0.05;
          end
          bAlphas = nan(B, 1);
          bBetas = nan(B, 1);
          bGammas = nan(B, 1);
          for i = 1:B
              bSample = datasample(sample, length(sample), 'replace', true);
              bParmHat = this.fitDist(bSample);
              bAlphas(i) = bParmHat(1);
              bBetas(i) = bParmHat(2);
              bGammas(i) = bParmHat(3);
          end
          
          % The index of the interval is chosen as in Efron and Tibshirani
          % (1993), p. 160.
          iLower = floor((B + 1) * (alpha / 2));
          iUpper = B + 1 - iLower;
          
          % Compute the estimators' standard deviations, see Eq. 2.3 in 
          % Efron and Tibshirani (1993).
          pStd = [std(bAlphas), std(bBetas), std(bGammas)];
          
          % Compute the 1 - alpha intervals based on bootstrap percentiles 
          % (see Efron and Tibshirani (1993), pp. 168). 
          sortedAlphas = sort(bAlphas);
          sortedBetas = sort(bBetas);
          sortedGammas = sort(bGammas);
          pCi = ...
              [sortedAlphas(iLower), sortedBetas(iLower), sortedGammas(iLower);
               sortedAlphas(iUpper), sortedBetas(iUpper), sortedGammas(iUpper)];
          
          parmHat = this.fitDist(sample); % Calling fitDist also sets the class' parameters.
          this.BootstrapParm = [bAlphas, bBetas, bGammas];
          this.ParameterSE = pStd;
          this.ParameterCI = pCi;
      end
      
      function f = pdf(this, x, alpha, beta, gamma)
          pdf = @(x, alpha, beta, gamma) (x > gamma) * ... % ensures that f(x<gamma) = 0.
              beta ./ alpha .* ((x - gamma) ./ ...
              alpha).^(beta - 1) .* exp(-1 .* ((x - gamma) ./ alpha).^beta);
      
          if nargin < 3
              f = pdf(x, this.Alpha, this.Beta, this.Gamma);
          else
              f = pdf(x, alpha, beta, gamma);
          end
      end
      
      function F = cdf(this, x)
          % Cumulative distribution function.
          F = (x > this.Gamma) * ... % ensures that F(x<gamma) = 0.
              (1 - exp(-1 .* ((x - this.Gamma) / this.Alpha).^this.Beta));
      end
      
      function x = icdf(this, p)
          % Inverse cumulative distribution function.
          x = wblinv(p, this.Alpha, this.Beta) + this.Gamma;
      end
      
      function x = drawSample(this, n)
          if n < 2
              n = 1;
          end
          p = rand(n, 1);
          x = this.icdf(p);
      end
      
      function val = negativeloglikelihood(this, x)
          % Negative log-likelihood value (as a metric of goodness of fit).
          val = sum(-log(pdf(x, this.Alpha, this.Beta, this.Gamma)));
      end
      
      function mae = meanabsoluteerror(this, sample, pi)
          % Mean absolute error (as a measure of goodness of fit).
          n = length(sample);
          if nargin < 3
              i = [1:n]';
              pi = (i - 0.5) ./ n;
          end
          xi = sort(sample);
          xhati = this.icdf(pi); % The prediction.
          mae = sum(abs(xi - xhati)) / n;
      end
      
      function ax = qqplot(this, sample, qqFig, qqAx, lineColor)
          if nargin > 2
              set(0, 'currentfigure', qqFig);
              set(qqFig, 'currentaxes', qqAx);
          else
              qqFig = figure();
          end
          if nargin < 4
              lineColor = [0 0 0];
          end
          n = length(sample);
          i = [1:n]';
          pi = (i - 0.5) ./ n;
          xi = sort(sample);
          xhati = this.icdf(pi); % The prediction.
          hold on
          plot(xhati, xi, 'kx'); 
          plot(xhati, xhati, 'color', lineColor, 'linewidth', 1.5);
          xlabel('Theoretical quantiles');
          ylabel('Ordered values');
      end
   end
end
