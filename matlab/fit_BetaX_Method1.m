function [Goodness,Paras,ci_1,ci_2]=fit_BetaX_Method1(alpha,gamma,x,y)
    % Fitting beta(x) curve using Method1, i.e., only estimates omega.
    % Input: alpha, intrinsic spreading power
    %        gamma, proportion of common neighbors
    %        x, number of exposures 
    %        y, retweeting probability
    % Output: Goodness, goodness-of-fit measures (rmse, rsquare, adjrsquare, sse, dfe)
    %         Paras, estimated parameters (omega)
    %         ci_1, 95% confidence intervals for the coefficient estimates
    %         ci_2, 95% prediction intervals for a new Y value at the specified X value
    [xData,yData]=prepareCurveData(x,y);
    fitfunc=[num2str(alpha),'*x.*((1-',num2str(gamma),').^(x.^omega))'];
    ft=fittype(fitfunc,'independent','x','dependent','y');
    opts=fitoptions('Method','NonlinearLeastSquares');
    opts.Algorithm='Trust-Region';
    opts.Lower=0;
    opts.Upper=Inf;
    opts.StartPoint=1.0;
    [fitresult,gof]=fit(xData,yData,ft,opts);
    ci_1=confint(fitresult,0.95);
    ci_2=predint(fitresult,x,0.95,'functional','on');
    Goodness=[gof.rmse,gof.rsquare,gof.adjrsquare,gof.sse,gof.dfe];
    Paras=fitresult.omega;
end