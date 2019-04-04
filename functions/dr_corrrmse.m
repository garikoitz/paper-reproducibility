function [N, Nold, Xm,Ym,rho,pval,rhom,pvalm,rmse,rmsem,rrmse,rrmsem,sdX,sdY,sdXm,sdYm,icc,iccm,CoV,CoVm] = dr_corrrmse(X,Y)
%corrrmse Some summary stats
%   Detailed explanation goes here
    oldmat            = [X, Y];
    Nold              = size(oldmat,1);
    maintainThisSubjs = (~sum(isnan(oldmat),2) ~=0);
    X                 = X(maintainThisSubjs, :);
    Y                 = Y(maintainThisSubjs, :);
    N                 = size(X,1);
    
    if Nold ~= N
        disp(['    WARNING: THERE WERE NANs, ' num2str(Nold-N) ' subjects from ' num2str(Nold) ' where deleted'])
    end
    
    if N == 0 
        Xm           = NaN;
        Ym           = NaN;
        rho          = NaN;
        pval         = NaN;
        rhom         = NaN;
        pvalm        = NaN;
        rmse         = NaN;
        rmsem        = NaN;
        rrmse        = NaN;
        rrmsem       = NaN;
        sdX          = NaN;
        sdY          = NaN;
        sdXm         = NaN;
        sdYm         = NaN;
        newmat       = NaN;
        keep         = NaN;
        icc          = NaN;
        iccm         = NaN;
        CoV          = NaN;
        CoVm         = NaN;
    else
        Xm           = mean(X,2,'omitnan');
        Ym           = mean(Y,2,'omitnan');
        [rho,pval]   = corr(X(:),Y(:),'type','Pearson','rows','pairwise','tail','both');
        [rhom,pvalm] = corr(Xm,Ym,'type','Pearson','rows','pairwise','tail','both');
        rmse         = sqrt(    mean(  (X(:)-Y(:)).^2, 'omitnan')    );
        rmsem        = sqrt(mean((Xm-Ym).^2, 'omitnan'));
        rrmse        = 100*rmse/((mean(Xm,'omitnan')+mean(Ym,'omitnan'))/2);
        rrmsem       = 100*rmsem/((mean(Xm,'omitnan')+mean(Ym,'omitnan'))/2);
        sdX          = std(X(:),'omitnan');
        sdY          = std(Y(:),'omitnan');
        sdXm         = std(Xm,'omitnan');
        sdYm         = std(Ym,'omitnan');
        CoV          = 100*( (sqrt(sum((X(:)-Y(:)).^2)/(2*length(X(:))))) / mean([X(:);Y(:)],'omitnan') );
        CoVm         = 100*( (sqrt(sum((Xm-Ym).^2)/(2*length(Xm)))) / mean([Xm;Ym]) );
        icc          = ICC([X(:),Y(:)],'3','1');
        iccm         = ICC([Xm,Ym],'3','1');
    end
end

