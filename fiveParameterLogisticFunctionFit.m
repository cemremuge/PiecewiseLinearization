function yf = fiveParameterLogisticFunctionFit(metricScore, mosScore)
%fiveParameterLogisticFunctionFit
% performs 5-parameter logistic function proposed by VQEG [1].
% 
% fiveParameterLogisticFunctionFit(metricScore, mosScore)
%   takes two input parameters as metricScore and mosScore, and uses these
%   values to assign an initial value to the parameters of the 5-parameter
%   logistic function. Then uses nonlinear regression functions (nlinfit and 
%   nlpredci) to optimize beta parameters.
%
% 
% [1] Video Quality Experts Group (VQEG)- "Final report on the validation 
% of objective models of video quality assessment",2003.
% 
% See also nlinfit, nlpredci.
% 
% C.Müge Bilsay
% 09.12.2024

x = metricScore;
y = mosScore;

% parameter initializaiton
b1 = max(y);
b2 = min(y);
b3 = mean(x);
b4 = 1;
b5 = 0;

% define the model
betaStartLF = [b1, b2, b3, b4, b5];
modelFiveParameterLF = @(b,x) b(1) .* (1/2-1 ./ (1+exp(b(2) .* (x - b(3))))) + b(4) .* x + b(5);

% find the initial beta values using nlinfit
modelFiveParameterLF = @(b,x) modelFiveParameterLF(b,x);
[bf, rw, jw] = nlinfit(x, y, modelFiveParameterLF, betaStartLF);

% confidence interval
xgrid = linspace(0, max(x), length(x));
[yf, delta] = nlpredci(modelFiveParameterLF, xgrid, bf, rw,jw);


end