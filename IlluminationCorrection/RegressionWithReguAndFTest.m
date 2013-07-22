function coeffs = RegressionWithReguAndFTest( InputInd, Vals )
%Setup parameters
FtestThresh = 0.9;
param.numThreads=-1;% all cores (-1 by default)
param.verbose=false; % verbosity, false by default
param.it0=10;       % frequency for duality gap computations
param.max_it=15000; % maximum number of iterations
param.lambda=0.1;   % regularization parameter
param.tol=1e-3;
param.intercept=true;
param.pos=false;
param.compute_gram=true;
% fprintf('\nFISTA + elstic net regularization\n');
param.loss='square';
param.regul='elastic-net';
param.lambda2=0.1;

% param.ista=false;
% param.subgrad=true;
% param.a=0.1;
% param.b=1000; % arbitrary parameters

W0=zeros(size(InputInd,2),size(Vals,2));

dontStop = 1;

RSSold = 0;
coeffs = W0;
%Loop to select the ideal number of features
while(dontStop)
    [intermCoeffs,optim_info]=mexFistaFlat(Vals,InputInd,W0,param);
    %intermCoeffs'
    if dontStop && max(max(optim_info))> 1e7
        disp('The optimization problem did not converge after 15k iters')
        dontStop = 0;
    end
    if dontStop
        if mod(dontStop,2)
            param.lambda  = param.lambda/2;
            param.lambda2 = param.lambda2/sqrt(2);
%             param.a=param.a/2;
%             param.b=param.b/2;
        end
        dontStop = dontStop + 1;
        %Compute the RSS
        RSS = Vals-InputInd*intermCoeffs;
        RSS = RSS'*RSS;
        %Skip f-test till one par then till the number of parameters goes up
        if sum(intermCoeffs~=0)==sum(coeffs~=0) || sum(coeffs~=0)==0
            coeffs = intermCoeffs;
            RSSold = RSS;
            continue
        end
        if RSSold<RSS
            dontStop = 0;
            continue
        else
            %Do F-test
            Fx = (RSSold-RSS)*(length(InputInd)-sum(intermCoeffs~=0)) / ...
                 ((sum(intermCoeffs~=0)-sum(coeffs~=0))*RSS);
            FCDF = fcdf1( Fx, (sum(intermCoeffs~=0)-sum(coeffs~=0)), ...
                             (length(InputInd)-sum(intermCoeffs~=0))      );
        end
        if FCDF<FtestThresh
            dontStop = 0;
            continue
        end
        RSSold = RSS;
        coeffs = intermCoeffs;
    end
    if sum(intermCoeffs==0)==0
%         fprintf('The solution is equivalent to L2 with no regularization\n');
        dontStop = 0;
    end
end
end