classdef TranslatedWeibull < handle
% The translated Weibull distribution, sometimes also simply called 
% "3-parameter Weibull distribution" is a parametric probability
% distribution.
%
% We use the parameterization and variables names that are also used in 
% doi.org/10.1016/j.coastaleng.2017.03.002 .
   properties
      Alpha
      Beta
      Gamma
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
          % Estimate the parameters of the distribution using maximum 
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
          F =  1 - exp(-1 .* ((x - this.Gamma) / this.Alpha).^this.Beta);
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
