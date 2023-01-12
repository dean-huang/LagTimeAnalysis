%% Complete analysis of E. coli

disp('refreshed')
clearvars
close all

% Change directory to main workspace.
dir_main = [fileparts(which('Ecoli_FullAnalysis.m')),'/'];
cd(dir_main);


%% Only need to run once:

%%% Filtering. Need to go into individual data folders.
% This is run separately for each data set. The resulting files will be
% stored in the exponential phase folder.

% Uses the median filter method: Takes the data, applies a median
% filter (200 kb), sorts the points based on how far they are from the 
% median filtered result, and then removes the most extreme 10%.

% Also filters out known problematic regions (due to sequencing-related
% biases) and divides the exponential phase data by the stationary phase
% data.

%{
mediums = {'LB','MM'};

for mm = 1:2
    
    medium = mediums{mm};
    % Directory of the bioproject:
    dir_bioproj = 'Data/Ecoli/';

    if isequal(medium,'LB')
        % Last four numbers of the ERR file name:
        exp_code = '3099';
        stat_code = '3100'; % Comment out if unused.
    elseif isequal(medium,'MM')
        exp_code = '3120';
        stat_code = '3121'; % Comment out if unused.
    end

    % Prefix of the file directories for exponential phase and stat phase:
    exp_dir = [medium,'_Exp'];
    stat_dir = [medium,'_Stat'];

    % Options below are 'oriWT' or 'oriZ'.
    ori = 'oriWT';

    % Options for filtering.
    filtwindow = 200; % Sets the size of the median filter.
    frac2remove = 0.10; % Sets the fraction of most extreme points to remove.


    %%% Loads the right files.

    cd([dir_main,dir_bioproj]);

    % Automates directory name and pos file name
    exp_posname = ['pos_',exp_code];
    phases = {'Exp'};
    if exist('stat_code','var')
        stat_posname = ['pos_',stat_code];
        phases = {'Exp','Stat'};
    end
    cd(exp_dir);

    L = 4641652;
    % nbin for 1kbp windows. Note that this is for the total
    % chromosome, not the left right legs individually.
    nbin = 4641;

    for kk = 1:numel(phases)
        phase = phases{kk};
        disp(phase)
        if isequal(phase,'Exp')
            data.(medium).Exp = nonzeros(load([exp_posname,'.mat']).pos);
        end
        if isequal(phase,'Stat')
            data.(medium).Stat = nonzeros(load(['../',stat_dir,'/',stat_posname,'.mat']).pos);
        end
    end

    if ~exist('Figs/','file')
        mkdir('Figs/')
    end

    % The name of the experiment consists of the medium, the chromosome type,
    % and the origin type.
    exptname = [medium,'_',ori]


    %%% Makes the 1 kb counts file. This will be filtered in the next section.

    disp(medium)
    for jj = 1:numel(phases)
        phase = phases{jj};
        disp(phase)

        tmp = data.(medium).(phase);
        tmp_1 = histcounts(tmp,1:1000:L)';
        counts_1 = tmp_1;
        save([exptname,'_',phase,'_1kb.mat'],'counts_1')

        % Plot unfiltered.
        h = figure();
        hold on;
        title([medium,' ',ori,' ',phase,' unfilt']);
        plot(counts_1);
        savefig(h,['Figs/',medium,'_',ori,'_',phase,'_unfilt']);

    end


    %%% Throw out the most extreme subset of points using the median filter.
    % Fraction is set by the frac2remove variable.

    disp(medium)
    for jj = 1:numel(phases)
        phase = phases{jj};
        disp(phase)
        unfilt = load([exptname,'_',phase,'_1kb.mat']).counts_1;

        % Do the filtering.
        counts_1 = medcutfilt(unfilt,filtwindow,frac2remove);

        % Plot filtered.
        h = figure();
        hold on;
        title([medium,' ',ori,' ',phase,' filt']);
        plot(counts_1);
        savefig(h,['Figs/',medium,'_',ori,'_',phase,'_filt_1kb_nohand']);

        % Compare with no filtering.
        h = figure();
        hold on;
        title([medium,' ',ori,' ',phase,' filter comparison']);
        plot(unfilt,'.','MarkerSize',3);
        plot(counts_1,'.','MarkerSize',3);
        savefig(h,['Figs/',medium,'_',ori,'_',phase,'_comp_filt']);

        % Save the data.
        save([exptname,'_',phase,'_filt_1kb_nohand.mat'],'counts_1');
    end

    % Normalize wrt stationary and filter the normalized data.
    if exist('stat_code','var')
        tmp_Exp = load([exptname,'_Exp_1kb.mat']).counts_1;
        tmp_Stat = load([exptname,'_Stat_1kb.mat']).counts_1;

        unfilt = tmp_Exp./tmp_Stat;

        counts_1 = medcutfilt(unfilt,filtwindow,frac2remove);

        % Plot filtered norm.
        h = figure();
        hold on;
        title([medium,' ',ori,' Norm filt']);
        plot(counts_1);
        savefig(h,['Figs/',medium,'_',ori,'_Norm_filt_1kb_nohand']);

        % Compare with no filtering.
        h = figure();
        hold on;
        title([medium,' ',ori,' ','Norm',' filter comparison']);
        plot(unfilt,'.','MarkerSize',3);
        plot(counts_1,'.','MarkerSize',3);
        savefig(h,['Figs/',medium,'_',ori,'_','Norm','_comp_filt']);

        save([exptname,'_','Norm','_filt_1kb_nohand.mat'],'counts_1');
    end
end

%}


%% NOTE: Everything below is an analysis for a single set of data.
% If you want to analyze the other sets, copy and paste from the Vibrio
% full analysis script and edit the code from there.

%% Also run once only:

disp('refreshed')
clearvars
close all

dir_main = [fileparts(which('Ecoli_FullAnalysis.m')),'/'];
cd(dir_main);

%%% Loads the right files into data structure and saves it.
%{
% Initialize experimental conditions here.
% These are all data from bioproject 2018_PRJEB25595, the 2018 Rudolph
% paper.

mediums = {'LB','MM'};

for mm = 1:2
    
    medium = mediums{mm};
    
    dir_bioproj = 'Data/Ecoli/';

    exp_dir = [medium,'_Exp'];

    if isequal(medium,'LB')
        exp_code = '3099';
    elseif isequal(medium,'MM')
        exp_code = '3120';
    end
    cd([dir_main,dir_bioproj,exp_dir]);
    
    ori = 'oriWT';
    handfilt = 'nohand';
    exptname = [medium,'_',ori];
    dataname = [exptname,'_Norm_filt_1kb_',handfilt,'.mat'];
    
    counts = load(dataname).counts_1(:,1);
    normcounts = counts/max(counts);
    data.(medium).(ori).counts = normcounts;
end

% Growth rates:
LBGT = 19.3;
data.LB.doubtimemins = LBGT;
data.LB.kG = log(2)/(LBGT*60);
data.LB.kGmins = log(2)/(LBGT);
MMGT = 68.8;
data.MM.doubtimemins = MMGT;
data.MM.kG = log(2)/(MMGT*60);
data.MM.kGmins = log(2)/(MMGT);

cd([dir_main,'Ecoli_FullAnalysis/']);
save(['data_',handfilt,'.mat'],'data')
%}


%%% Shift data so that oriC is in the middle. Also makes 10kb bins and std.
%{

plotfigs = 1;
do_rDNA = 0;

ori = 'oriWT';
handfilt = 'nohand';

data = load(['data_',handfilt,'.mat']).data;

% Secondary binning of the 1kbp data.
bin_var = 10; % in kbp.
data.bin_var = bin_var;

L_full = 4641652;
data.L_full = L_full;

oriC_pos = [3925744, 3925975];
oriC = floor(mean(oriC_pos)/1000);
data.oriC1k = oriC;
oriCfull = floor(mean(oriC_pos));
data.oriC = oriCfull;

% Anonymous function for shifting the right arm by oriC location
shift = @(x) mod(x-oriCfull,L_full);

% rDNA locations right arm
rrs_r_label = {'rrsC';'rrsA';'rrsB';'rrsE';'rrsH'};
rrs_r_pos = [ ...
    3941808, 3943349; ... % rrsC
    4035531, 4037072; ... % rrsA
    4166659, 4168200; ... % rrsB
    4208147, 4209688; ... % rrsE
    223771, 225312];      % rrsH
rrs_r_mean = mean(rrs_r_pos,2);
rrs_r = shift(rrs_r_mean);
% rDNA locations left arm
rrs_l_label = {'rrsD';'rrsG'};
rrs_l_pos = [ ...
    3427221, 3428762; ... % rrsD
    2729616, 2731157];    % rrsG
rrs_l_mean = mean(rrs_l_pos,2);
rrs_l = rrs_l_mean-oriCfull; % Left arm shifted just by subtracting by oriC.

data.rrs = [rrs_l;rrs_r];
data.rrs_label = [rrs_l_label;rrs_r_label];

mediums = {'LB','MM'};

for mm = 1:2
    
    medium = mediums{mm};
    counts = data.(medium).(ori).counts;
    logy = log(counts);
    
    sigma_counts = Get_SigmaEst(counts);
    sigma_logy = Get_SigmaEst(logy);
    
    oriat0 = circshift(counts,-oriC);
    
    % The following changes the x values such that the oriC is at 0 and
    % the left leg has negative values. The 1kbp bins are placed at half
    % kbp locations. i.e. what was originally bin 1 goes to 0.5, bin 2 goes
    % to 1.5, and bin (half+1) goes to -0.5, etc.
    x_prior = (1:numel(oriat0))';
    xhalf = ceil(numel(x_prior)/2);
    x_new = x_prior-xhalf-0.5;
    counts_shifted = circshift(oriat0,xhalf);
    nnanind = ~isnan(counts_shifted);
    y = counts_shifted(nnanind);
    x = x_new(nnanind);
    
    data.(medium).(ori).y = y;
    data.(medium).(ori).x = x*1000;
    data.(medium).(ori).sigma_counts = sigma_counts;
    data.(medium).(ori).sigma_logy = sigma_logy;
    
    % Makes a binned version using binsize in kbp.
    x_var_edges = floor(min(x)/10)*10:10:ceil(max(x)/10)*10;
    x_var = mean([x_var_edges(1:end-1);x_var_edges(2:end)],1);
    y_var = zeros(size(x_var));
    y_var_sigma = zeros(size(x_var));
    for jj = 1:numel(x_var_edges)-1
        left = x_var_edges(jj);
        right = x_var_edges(jj+1);
        y_var(jj) = mean(counts_shifted(x_new >= left & x_new < right),'omitnan');
        % The error for each of these binned points is just the standard
        % deviation of the individual 1kbp bin points.
        y_var_sigma(jj) = std(counts_shifted(x_new >= left & x_new < right),'omitnan');
    end
    nnanind_var = ~isnan(y_var);
    x_var = x_var(nnanind_var)*1000;
    y_var = y_var(nnanind_var);
    y_var_sigma = y_var_sigma(nnanind_var);
    
    data.(medium).(ori).y_var = y_var;
    data.(medium).(ori).x_var = x_var;
    data.(medium).(ori).y_var_sigma = y_var_sigma;
    
    if plotfigs == 1
        figure()
        title([medium,' ',ori])
        hold on;
        plot(x,y)
        set(gca,'YScale','log')
        if do_rDNA == 1
            for ii = 1:numel(data.rrs)
                rlines = xline(data.rrs(ii)/1000,'--r',data.rrs_label(ii),...
                    'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
            end
        end
        figure()
        title([medium,' ',ori,' 10 kb'])
        hold on;
        plot(x_var,y_var)
        set(gca,'YScale','log')
    end
end

% Save the data with shifted versions.
save('data_shifted.mat','data')
%}


%%% Get lag time data from single slope fit. Saved into data structure.
%{

debug = 1;
do_legfixplot = 1;
do_shade = 1;
do_normbyfit = 1;

% Load data.
data = load('data_shifted.mat').data;

mediums = {'LB','MM'};

for mm = 1:2
    
    medium = mediums{mm};

    % Get kG.
    kGmins = data.(medium).kGmins;
    
    % Get relevant data.
    x = data.(medium).oriWT.x';
    y = data.(medium).oriWT.y';
    logy = log(y);
    sigma = data.(medium).oriWT.sigma_logy;
    
    % Get control point locations.
    xpts = [min(x),0,max(x)];
    guess = [1,1,1,xpts];
    
    % Do fit on log data.
    [ypts,~,logyfit] = Fit_PiecewiseLinear_nLeg_SetX(x,logy,guess,kGmins,sigma,debug);
    
    % Normalize with respect to fit.
    if do_normbyfit == 1
        ymax = max(logyfit);
        logyfit = logyfit-ymax;
        logy = logy-ymax;
        y = y./exp(ymax);
    end
    
    % Plot.
    figure();
    plot(x,logy);
    hold on;
    plot(x,logyfit,'LineWidth',3);
    title([medium,' oriWT']);
    
    % Get lag time from fit data.
    tau = -(kGmins^(-1))*logyfit;
    data.(medium).oriWT.tau = tau';

end

% Save the data with the lag times.
save('data_shifted.mat','data');
%}


%%% Do AIC.
% Chooses 32 legs for full LB data, but 16 legs is pretty close.
%{

do_plt = 0;
do_sigmaLL = 1;
debug = 0;

% Number of hypotheses. Highest number of legs is 2^(ncond). Lowest is 2.
n_hyp = 7;

data = load('data_shifted.mat').data;

mediums = {'LB','MM'};

for mm = 1:2
    
    medium = mediums{mm};
    
    logy = log(data.(medium).oriWT.y);
    x = data.(medium).oriWT.x;
    kG = data.(medium).kG;
    sigma = data.(medium).oriWT.sigma_logy;
    
    % sigmaLL uses the fact that we care about fork vel as our
    % observable, so we would have to double count the variance when
    % dealing with the log likelihood estimate.
    if do_sigmaLL == 1
        sigmaLL = sqrt(2)*sigma;
    else
        sigmaLL = sigma;
    end
    
    LL = zeros(n_hyp,1);
    
    % +1 because the number of fit parameters is number of legs +1.
    nparam_max = 2^(n_hyp)+1;
    params = nan(n_hyp,nparam_max);
    param_errs = nan(n_hyp,nparam_max);
    
    for nn = 1:n_hyp
        % Number of segments in hypothesis. nn=1 and nleg=2 for null hypothesis.
        n_leg = 2^nn;
        legsize = min(2*[floor((max(x)-0)/n_leg),floor((0-min(x))/n_leg)]);
        xpts = [fliplr(0:-legsize:floor(min(x))),legsize:legsize:ceil(max(x))];
        ypts = zeros(1,numel(xpts));
        guess = [ypts,xpts];
        
        % Get rid of the points that lie outside the xpts domain.
        outside = x < min(xpts) | x > max(xpts);
        logy(outside) = [];
        x(outside) = [];
        
        % Continuous piecewise fit.
        [param,param_err,yfit] = Fit_PiecewiseLinear_nLeg_SetX(x,logy,guess,kG,sigma,debug);
        params(nn,1:n_leg+1) = param;
        param_errs(nn,1:n_leg+1) = param_err;
        
        if do_plt
            figure();
            plot(x,logy,'b.');
            hold on;
            plot(x,yfit,'-','LineWidth',3);
            title([medium,' ',num2str(n_leg),' legs']);
            
            errorbar(xpts,param,param_err,'o','MarkerSize',10,'LineWidth',3);
            
            pause(1)
            
        end
        
        LL(nn) = - sum((logy-yfit).^2)/(2*sigmaLL^2);
        
    end
    disp('log likelihood');
    disp(LL);
    
    [logLambda, CDF_chisq, Pval, nparams, AIC] = deal(nan(n_hyp,1));
    for nn = 2:n_hyp
        nleg = 2^nn;
        chidim = (nleg+1)-(nleg/2+1); % Wilks' theorem
        % Log of the likelihood ratio is equal (up to ambiguous factors of 2)
        % to the chi squared distribution of dimension = dim(alt. hyp.) - dim(null hyp.)
        % Number of parameters is initial intercept + slope of each leg
        % Null hypothesis has intercept + slope = dim 2
        % Two leg hypothesis has intercept + slope + slope = dim 3
        
        logLambda(nn) = LL(nn)-LL(nn-1); % each subsequent condition is compared to the one prior
        % (2 legs vs. 1 leg, 4 legs vs. 2 legs, 8 legs vs. 4 legs, etc.)
        CDF_chisq(nn) = chi2cdf(2*logLambda(nn),chidim);
        Pval(nn) = 1-CDF_chisq(nn);
        
        % AIC
        nparams(nn) = nleg+1;
        AIC(nn) = nparams(nn) - LL(nn);
        
    end
    AIC(1) = 3 - LL(1);
    
    figure();
    plot(AIC-min(AIC));
    title([medium,' Delta AIC']);
    legs = 2.^(1:n_hyp);
    xticklabels(strsplit(num2str(2.^(1:n_hyp))));
    
    % disp('CDF');
    % disp(CDF_chisq);
    disp('Pval');
    disp(Pval);
    disp('AIC');
    disp(AIC);
    
    % Use AIC best choice to plot and with control point errors
    
    [~,nhyp_best] = min(AIC);
    nleg_best = 2^(nhyp_best)
    
    legsize = min(2*[floor((max(x)-0)/nleg_best),floor((0-min(x))/nleg_best)]);
    xpts = [fliplr(0:-legsize:floor(min(x))),legsize:legsize:ceil(max(x))];
    guess = [ones(1,numel(xpts)),xpts];
    outside = x < min(xpts) | x > max(xpts);
    logy(outside) = [];
    x(outside) = [];
    [param,param_err,yfit,rep,rep_err] = Fit_PiecewiseLinear_nLeg_SetX(x,logy,guess,kG,sigma,debug);
    
    figure();
    plot(x,logy,'b.');
    hold on;
    plot(x,yfit,'-','LineWidth',3);
    errorbar(xpts,param,param_err,'o','MarkerSize',10,'LineWidth',3);
    plot(x,yfit,'-','LineWidth',3);
    title([medium, ' fit'])
    figure();
    errorbar(xpts(2:end),rep,rep_err,'o','MarkerSize',10,'LineWidth',3);
    title([medium, ' rep'])
end
%}


%% Figures:

%%% Make lag time figure. Currently uses 10 kb bins.
%{
debug = 0;
do_shade = 1;
do_plot1kbp = 0;
do_lagtime = 1;
do_fixsigs = 0;
do_singleslope = 0;
do_legfixplot = 0;
do_normbyfit = 1;

if do_fixsigs == 1
    do_singleslope = 0;
end

data = load('data_shifted.mat').data;

mediums = {'LB','MM'};

for mm = 1:2
    
    medium = mediums{mm};
    
    h = figure(100+mm);
    hold on;
    box on;
    
    kGmins = data.(medium).kGmins;
    kG = data.(medium).kG;
    
    % The following is for plotting 1kbp points. If you want to plot other
    % sizes, set it to 0 and modify the second portion of the loop.
    if do_plot1kbp == 1
        binsize = 1000;
        x1 = data.(medium).oriWT.x';
        y1 = data.(medium).oriWT.y';
        sig1 = data.(medium).oriWT.sigma_counts;
        %ymax = max(y1);
        %y1 = y1./ymax;
        logy1 = log(y1);
    else
        binsize = data.bin_var*1000;
        x1 = data.(medium).oriWT.x_var;
        y1 = data.(medium).oriWT.y_var;
        sig1 = data.(medium).oriWT.y_var_sigma;
        %ymax = max(y1);
        %y1 = y1./ymax;
        logy1 = log(y1);
    end
    
    % To make the fit be a single slope for all data, including both chrs.
    if do_singleslope == 1
        
        disp('NEED TO REWRITE SINGLE SLOPE CODE')
        
    else
        Ndata1 = numel(y1);
        xpts1 = [min(x1),0,max(x1)];
        guess1 = [1,1,1,xpts1];
        [ypts1,~,yfit1] = Fit_PiecewiseLinear_nLeg_SetX(x1,log(y1),guess1,[],[],debug);
        
        % To fix the standard deviations being 0 for the 10 kbp bins that had only
        % a single 1 kbp point within, let the standard deviation be equal to the
        % difference between the 10 kbp bin value and the single slope fit value
        % (which is effectively the mean).
        if do_fixsigs == 1
            for ii = 1:numel(sig1)
                if sig1(ii) == 0
                    disp(ii)
                    sig1(ii) = abs(y1(ii)-exp(yfit1(ii)));
                    disp(sig1(ii))
                end
            end
            data.(medium).oriWT.y_var_sigma = sig1;
            save('data_shifted.mat','data')
        end
        
        % Normalize by the yfit so that lag time of 0 is set by the fit rather
        % than the max of y.
        if do_normbyfit == 1
            logymax = max(ypts1);
            yfit1 = yfit1-logymax;
            ypts1 = ypts1-logymax;
            logy1 = logy1-logymax;
            y1 = y1./exp(logymax);
        end
        
        L3 = plot(xpts1,exp(ypts1),'-','Color',[.929,.694,.125],'LineWidth',1.5);
    end
    
    figure(h);
    L1 = semilogy(x1,y1,'-','Color',[0,.447,.741]);
    
    if do_shade == 1
        patch([x1,fliplr(x1)],[y1-sig1,fliplr(y1+sig1)],...
            [0,.447,.741],'FaceAlpha',0.15,'EdgeColor','none');
    end
    
    %title('LB oriWT')
    xlabel('Position relative to{\it oriC}: L (bp)')
    ylabel('Locus copy number: N_L/N_0')
    legendnames = {'10 kb sampling','Fit'};
    legend([L1,L3],legendnames,'Location','northeast')
    if mm == 1
        YLimL = [0.3,1.1];
    elseif mm == 2
        YLimL = [0.5,1.1];
    end
    ylim(YLimL)
    set(gca,'YScale','log')
    
    if do_lagtime == 1
        yyaxis right
        ylabel('Lag time: \tau_L (min)');
        set(get(gca,'YLabel'),'Rotation',-90,'VerticalAlignment','Bottom')
        %set(gca,'YLim',YLimR)
        %set(gca,'YTickLabel',YTickLabelR)
        YLimR = -kGmins^-1*log(YLimL(end:-1:1));
        ylim(YLimR)
        set( gca, 'ydir', 'reverse' );
    end
    
    doPageFormat([5,3])
    print([medium,'_oriWT_CopyNumLagTime.pdf'],'-dpdf')
    savefig(h,[medium,'_oriWT_CopyNumLagTime'])

end

if do_fixsigs == 1
    disp('Please change do_fixsigs to 0 and run this again.')
end
%}


%%% Make fork velocity figure for single data set.
% Spacing set by xdel. Can also use number of legs to set.
%{
debug = 0;
do_legend = 1;
do_arrows = 1;
do_byxdel = 1; % 1 is with a set delta x (e.g. 200 kb), 0 is by number of legs (e.g. 32)
do_rDNA = 0;
do_rDNAlabel = 0;

data = load('data_shifted.mat').data;

medium = 'LB';

if do_byxdel == 1
    xdel = 2e5
else
    numlegs = 32
    x = data.(medium).oriWT.x;
    xdel = (max(x)-min(x))/(numlegs+1)
end
orioffset = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = figure();
hold on;
ax = gca;
hold on;
box on;

L = gobjects(4,1);
colors = {[0,.447,.741];[.85,.325,.098];[0,.447,.741];[.85,.325,.098]};

x = data.(medium).oriWT.x;
xsize = floor((max(x)-0)/xdel)+floor((0-min(x))/xdel);
reps = zeros([xsize,1]);

ii = 1;
jj = 0;

logy = log(data.(medium).oriWT.y);
x = data.(medium).oriWT.x;
kG = data.(medium).kG;
sigma = data.(medium).oriWT.sigma_logy;

for LR = 1:2 % Left or right leg index.
    
    if LR == 1
        xpts = fliplr(-orioffset:-xdel:min(x));
    else
        xpts = orioffset:xdel:max(x);
    end
    ypts = ones(size(xpts));
    
    guess = [ypts,xpts];
    
    [param,param_err,yfit,rep,rep_err] = Fit_PiecewiseLinear_nLeg_SetX(x,logy,guess,kG,sigma,debug);
    rep = abs(rep);
    rep_err = abs(rep_err);
    
    xmids = mean([xpts(1:end-1);xpts(2:end)],1);
    
    reps(jj+1:jj+numel(rep)) = rep;
    jj = jj+numel(rep);
    
    L(ii) = plot(xmids,rep,'.-','Color',colors{ii},'LineWidth',0.5,'MarkerSize',12);
    
    patch([xmids,fliplr(xmids)], [rep-rep_err,fliplr(rep+rep_err)], colors{ii}, 'FaceAlpha',0.3, 'EdgeColor','none');
end

if isequal(medium,'LB')
    maxylim = 2500;
else
    maxylim = 2500;
end

% Setting limits.
ylim([0,maxylim]);

% Making ytick labels be in kb/s.
yts = 0:500:maxylim;
ytls = yts/1000;
set(ax,'YTick',yts,'YTickLabel',ytls);

% Making xtick labels be in Mb.
xts = -3e6:1e6:3e6;
ax.XTick = xts;
xtls = strsplit(num2str(-3:1:3));
ax.XTickLabel = xtls;

% oriC line.
xline(0,'--','{\it oriC}','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');

% Mean LB rates.
repmean = yline(mean(reps),'--','Color',[0,.447,.741]);

% xlabel and ylabel.
xlabel(ax,'Position relative to{\it oriC}: L (Mb)');
ylabel(ax,'Fork velocity: v (kb/s)');

axpos = ax.Position;

if do_rDNA == 1
    for ii = 1:numel(data.rrs)
        if do_rDNAlabel == 1
            rlines = xline(data.rrs(ii),':r',data.rrs_label(ii),...
                'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
        else
            rlines = xline(data.rrs(ii),':r');
        end
    end
end

if do_legend == 1
    legendnames = {medium,'Mean'};
    lg = legend([L(1);repmean],legendnames,'Location','north');
    if do_rDNA == 1
        legendnames = {medium,'Mean','rDNA'};
        lg = legend([L(1);repmean;rlines(1)],legendnames,'Location','north');
    end
    lgpos = lg.Position;
    lg.Location = 'none';
    lgpos(2) = axpos(2)+axpos(4)-0.01;
    lg.Position = lgpos;
    lg.Orientation = 'horizontal';
    lg.FontSize = 10;
    legend('boxoff');
end

if do_arrows == 1
    xlimtot = get(gca,'XLim');
    xlimdom = xlimtot(2)-xlimtot(1);
    xdist = 0-xlimtot(1);
    yloc = axpos(2)+axpos(4);
    arry = yloc - 0.04 + [0,0];
    
    for LR = 1:2
        arrx = axpos(1)+axpos(3)*(xdist)/(xlimdom);
        arrx = arrx + (-1).^(LR)*(0.02 + [0,0.05]);
        annotation('arrow',arrx,arry);
        arbarx = [arrx(1),arrx(1)];
        arbary = arry(1)+[-0.01,0.01];
        annotation('line',arbarx,arbary);
    end
end

doPageFormat([5,3])
print([medium,'_forkvel.pdf'],'-dpdf')
savefig(h,[medium,'_forkvel'])
%}


%%% Make LBvMM fork velocity figure.
% Spacing set by xdel. Can also use number of legs to set.
%{
debug = 0;
do_legend = 1;
do_arrows = 1;
do_byxdel = 1; % 1 is with a set delta x (e.g. 200 kb), 0 is by number of legs (e.g. 32)
do_rDNA = 1;
do_rDNAlabel = 0;
do_diffplot = 0;

data = load('data_shifted.mat').data;

if do_byxdel == 1
    xdel = 2e5
else
    numlegs = 16
    LBx = data.LB.oriWT.x;
    xdel = (max(LBx)-min(LBx))/(numlegs+1)
end
orioffset = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = figure(104);
hold on;
ax = gca;
hold on;
box on;

L = gobjects(4,1);
colors = {[0,.447,.741];[.85,.325,.098];[0,.447,.741];[.85,.325,.098]};

LBx = data.LB.oriWT.x;
LBxsize = floor((max(LBx)-0)/xdel)+floor((0-min(LBx))/xdel);
LBreps = zeros([LBxsize,1]);
MMx = data.MM.oriWT.x;
MMxsize = floor((max(MMx)-0)/xdel)+floor((0-min(MMx))/xdel);
MMreps = zeros([MMxsize,1]);

mediums = {'LB','MM'};

ii = 0;
for mm = 1:2
    
    medium = mediums{mm};
    
    jj = 0;
    ii = ii+1;
    
    logy = log(data.(medium).oriWT.y);
    x = data.(medium).oriWT.x;
    kG = data.(medium).kG;
    sigma = data.(medium).oriWT.sigma_logy;
    
    for LR = 1:2 % Left or right leg index.
        
        if LR == 1
            xpts = fliplr(-orioffset:-xdel:min(x));
        else
            xpts = orioffset:xdel:max(x);
        end
        ypts = ones(size(xpts));
        
        guess = [ypts,xpts];
        
        [param,param_err,yfit,rep,rep_err] = Fit_PiecewiseLinear_nLeg_SetX(x,logy,guess,kG,sigma,debug);
        rep = abs(rep);
        rep_err = abs(rep_err);
        
        xmids = mean([xpts(1:end-1);xpts(2:end)],1);
        
        if isequal(medium,'LB')
            LBreps(jj+1:jj+numel(rep)) = rep;
        else
            MMreps(jj+1:jj+numel(rep)) = rep;
        end
        jj = jj+numel(rep);
        
        L(ii) = plot(xmids,rep,'.-','Color',colors{ii},'LineWidth',0.5,'MarkerSize',12);
        
        patch([xmids,fliplr(xmids)], [rep-rep_err,fliplr(rep+rep_err)], colors{ii}, 'FaceAlpha',0.3, 'EdgeColor','none');
    end
end

% Setting limits.
ylim([0,2500]);

% Making ytick labels be in kb/s.
yts = 0:500:2500;
ytls = yts/1000;
set(ax,'YTick',yts,'YTickLabel',ytls);

% Making xtick labels be in Mb.
xts = -3e6:1e6:3e6;
ax.XTick = xts;
xtls = strsplit(num2str(-3:1:3));
ax.XTickLabel = xtls;

% oriC line.
xline(0,'--');
text(0,ax.YLim(2),'{\itoriC1}','HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',9);

% Mean rates.
LBmean = yline(mean(LBreps),'--','Color',[0,.447,.741]);
MMmean = yline(mean(MMreps),'--','Color',[.85,.325,.098]);

% xlabel and ylabel.
xlabel(ax,'Position relative to{\it oriC}: L (Mb)');
ylabel(ax,'Fork velocity: v (kb/s)');

axpos = ax.Position;

if do_rDNA == 1
    for ii = 1:numel(data.rrs)
        if do_rDNAlabel == 1
            rlines = xline(data.rrs(ii),':k',data.rrs_label(ii),...
                'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
        else
            rlines = xline(data.rrs(ii),':k');
        end
    end
end

if do_legend == 1
    legendnames = {'LB','MM','LB Mean', 'MM Mean'};
    lg = legend([L(1:2);LBmean;MMmean],legendnames,'Location','north');
    if do_rDNA == 1
        legendnames = {'LB','MM','Mean','rDNA'};
        lg = legend([L(1:2);LBmean;rlines(1)],legendnames,'Location','north');
    end
    lgpos = lg.Position;
    lg.Location = 'none';
    lgpos(2) = axpos(2)+axpos(4)-0.02;
    lg.Position = lgpos;
    lg.Orientation = 'horizontal';
    lg.FontSize = 10;
    legend('boxoff');
end

if do_arrows == 1
    xlimtot = get(gca,'XLim');
    xlimdom = xlimtot(2)-xlimtot(1);
    xdist = 0-xlimtot(1);
    yloc = axpos(2)+axpos(4);
    arry = yloc - 0.04 + [0,0];
    
    for LR = 1:2
        arrx = axpos(1)+axpos(3)*(xdist)/(xlimdom);
        arrx = arrx + (-1).^(LR)*(0.02 + [0,0.05]);
        annotation('arrow',arrx,arry);
        arbarx = [arrx(1),arrx(1)];
        arbary = arry(1)+[-0.01,0.01];
        annotation('line',arbarx,arbary);
    end
end

doPageFormat([5,3])
print('LBvMM_forkvel.pdf','-dpdf')
savefig(h,'LBvMM_forkvel')

if do_diffplot == 1
    LB_diff = LBreps-mean(LBreps);
    MM_diff = MMreps-mean(MMreps);
    figure();
    plot(LB_diff);
    hold on;
    plot(MM_diff);
    legend({'LB diff','MM diff'});
end
%}


%%% Make velocity vs. lag time figure (VelLag).
% Only for one medium. You can choose whether the x interval differences
% are fixed or the lag time intervals.
%{
debug = 0;
do_legend = 1;
% Fixed xdels rather than taudels. The x interval differences are fixed
% rather than the lag time intervals. Remember to change them in the for
% loop below.
do_fixxdel = 1;
do_bynumlegs = 0; % 0 is with a set delta x (e.g. 200 kb), 1 is by number of legs (e.g. 32)

medium = 'MM';

data = load('data_shifted.mat').data;

if do_fixxdel == 1
    if do_bynumlegs == 1
        numlegs = 16
        xdel = (max(x)-min(x))/(numlegs+1)
    else
        xdel = 2e5
    end
    x = data.(medium).oriWT.x;
    xhalf = sum(x<0)+1;
    tau = data.(medium).oriWT.tau;
    xpts = 0:xdel:x(end);
    taupts = (tau(end)-tau(xhalf))*xpts/(x(end)-x(xhalf));
else
    taudel = 4
    tau1min = 0;
    tau1max = 44;
    x = data.(medium).oriWT.x;
    xhalf = sum(x<0)+1;
    tau = data.(medium).oriWT.tau;
    taupts = tau1min:taudel:max(tau);
    xpts = (x(end)-x(xhalf))*taupts/(tau(end)-tau(xhalf));
end
orioffset = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = figure(105);
hold on;
box on;

kG = data.(medium).kG;

L = gobjects(4,1);
colors = {[0,.447,.741];[.85,.325,.098];[0,0,1];[1,0,0]};

ii = 0;
x = data.(medium).oriWT.x;
logy = log(data.(medium).oriWT.y);
sigma = data.(medium).oriWT.sigma_logy;
tau = data.(medium).oriWT.tau;

xhalf = sum(x<0)+1;
xR = x(xhalf:end);

figure(1);
hold on;
plot(x,tau);
for ylines = 1:numel(taupts)
    yline(taupts(ylines),':');
end
figure(105);

taumids = mean([taupts(1:end-1);taupts(2:end)],1);

for LR = 1:2 % Left or right leg index.
    
    ii = ii+1;
    
    if LR == 1
        xptsLR = fliplr(-xpts);
        taumidsLR = fliplr(taumids);
    elseif LR == 2
        xptsLR = xpts;
        taumidsLR = taumids;
    end
    
    ypts = ones(size(xptsLR));
    
    guess = [ypts,xptsLR];
    
    [~,~,yfit,rep,rep_err] = Fit_PiecewiseLinear_nLeg_SetX(x,logy,guess,kG,sigma,debug);
    rep = abs(rep);
    rep_err = abs(rep_err);
    
    L(ii) = plot(taumidsLR,rep,'.-','Color',colors{ii},'LineWidth',1.5,'MarkerSize',12);
    
    patch([taumidsLR,fliplr(taumidsLR)], [rep-rep_err,fliplr(rep+rep_err)], colors{ii}, 'FaceAlpha',0.3, 'EdgeColor','none');
end

ylim([0,2500]);
if do_fixxdel == 0
    xlim([tau1min,tau1max]);
end

% Making ytick labels be in kb/s.
yts = 0:500:2500;
ytls = yts/1000;
set(gca,'YTick',yts,'YTickLabel',ytls);

% xlabel and ylabel.
xlabel('Lag time: \tau_L (min)');
ylabel('Fork velocity: v (kb/s)');

axpos = get(gca,'Position');

if do_legend == 1
    legendnames = {'L','R'};
    lg = legend(L(1:2),legendnames,'Location','northeast');
end

doPageFormat([5,3])
print([medium,'_VelLag.pdf'],'-dpdf')
savefig(h,[medium,'_VelLag'])
%}


%%% Make velocity vs. lag time figure to compare media (LBvMM_VelLag).
% For comparison. You can choose whether the x interval differences
% are fixed or the lag time intervals.
%{
debug = 0;
do_legend = 1;
% Fixed xdels rather than taudels. The x interval differences are fixed
% rather than the lag time intervals. Remember to change them in the for
% loop below.
do_fixxdel = 0;
do_bynumlegs = 0; % 0 is with a set delta x (e.g. 200 kb), 1 is by number of legs (e.g. 32)
do_diff = 0; % 1 makes the plot a difference from the mean.

data = load('data_shifted.mat').data;

mediums = {'LB','MM'};

ii = 0;
legendnames = cell(1,4);
L = gobjects(4,1);
%colors = {[0,.447,.741];[.85,.325,.098];[0,0,1];[1,0,0]};
colors = {rgb(255, 171, 118);rgb(255, 99, 99);rgb(107, 203, 119);rgb(77, 150, 255)};

for mm = 1:2
    
    medium = mediums{mm};
    
    if do_fixxdel == 1
        if do_bynumlegs == 1
            numlegs = 16
            xdel = (max(x)-min(x))/(numlegs+1)
        else
            xdel = 2e5
        end
        x = data.(medium).oriWT.x;
        xhalf = sum(x<0)+1;
        tau = data.(medium).oriWT.tau;
        xpts = 0:xdel:x(end);
        taupts = (tau(end)-tau(xhalf))*xpts/(x(end)-x(xhalf));
    else
        taudel = 4
        tau1min = 0;
        tau1max = 44;
        x = data.(medium).oriWT.x;
        xhalf = sum(x<0)+1;
        tau = data.(medium).oriWT.tau;
        taupts = tau1min:taudel:max(tau);
        xpts = (x(end)-x(xhalf))*taupts/(tau(end)-tau(xhalf));
    end
    orioffset = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    h = figure(105);
    hold on;
    box on;
    
    kG = data.(medium).kG;
    
    x = data.(medium).oriWT.x;
    logy = log(data.(medium).oriWT.y);
    sigma = data.(medium).oriWT.sigma_logy;
    tau = data.(medium).oriWT.tau;
    
    xhalf = sum(x<0)+1;
    xR = x(xhalf:end);
    
    figure(1);
    hold on;
    plot(x,tau);
    for ylines = 1:numel(taupts)
        yline(taupts(ylines),':');
    end
    figure(105);
    
    taumids = mean([taupts(1:end-1);taupts(2:end)],1);
    
    for LR = 1:2 % Left or right leg index.
        
        ii = ii+1;
        
        if LR == 1
            xptsLR = fliplr(-xpts);
            taumidsLR = fliplr(taumids);
        elseif LR == 2
            xptsLR = xpts;
            taumidsLR = taumids;
        end
        
        ypts = ones(size(xptsLR));
        
        guess = [ypts,xptsLR];
        
        [~,~,yfit,rep,rep_err] = Fit_PiecewiseLinear_nLeg_SetX(x,logy,guess,kG,sigma,debug);
        rep = abs(rep);
        rep_err = abs(rep_err);
        
        if do_diff == 1
            rep = rep-mean(rep);
            L(ii) = plot(taumidsLR,rep,'.-','LineWidth',1.5,'MarkerSize',12);
        else
            L(ii) = plot(taumidsLR,rep,'.-','Color',colors{ii},'LineWidth',1.5,'MarkerSize',12);
            patch([taumidsLR,fliplr(taumidsLR)], [rep-rep_err,fliplr(rep+rep_err)], colors{ii}, 'FaceAlpha',0.3, 'EdgeColor','none');
        end
        
        if LR == 1
            legendnames{ii} = [medium,' L'];
        else
            legendnames{ii} = [medium,' R'];
        end
    end
end
    
if do_fixxdel == 0
    xlim([tau1min,tau1max]);
end

ylim([0,2500]);
% Making ytick labels be in kb/s.
yts = 0:500:2500;
if do_diff == 1
    ylim([-500,500]);
    yts = -500:500:500;
    yline(0,':')
end
ytls = yts/1000;
set(gca,'YTick',yts,'YTickLabel',ytls);

% xlabel and ylabel.
xlabel('Lag time: \tau_L (min)');
ylabel('Fork velocity: v (kb/s)');

axpos = get(gca,'Position');

if do_legend == 1
    figure(105);
    lg = legend(L, legendnames, 'Location', 'northeast');
end

doPageFormat([5,3])
print(['LBvMM','_VelLag.pdf'],'-dpdf')
savefig(h,['LBvMM','_VelLag'])
%}


%% Make VelPosLag panel (super figure).
% Top left: logN+lag time vs. L. Full L.
% Bottom left: v vs. L, sharing x axis with top left. 
% Top right: Distance vs. time. Lag time starts at 0.
% Bottom right: v vs. time, sharing x axis with top right.
% crtS only on the time plot and extending through to both right hand plots.
% Color scheme: 1L is blue, 1R is red, 2L is cyan, 2R is magenta. 
% If not using Vc 2chr, then cyan and magenta are used for other media. 
% For shaded graphs, no lines are plotted.
%%{
debug = 0;
do_arrows = 1;
do_legend = 1;
do_dotdash = 0;
do_varbinsize = 0;
do_ymaxnorm = 0;
do_ABCD = 1;
do_ypts = 0;
do_yptserrbar = 0;

filename = 'LBMM_oriWT_VelPosLagPanel';

conditions = {'1L','1R','2L','2R'};
colors = {'b','r',rgb(0,102,168),rgb(174,0,127)};
if do_dotdash == 1
    marker = '.-';
else
    marker = '-';
end

% Start plotting.
h = figure(100);
%clf;
p = panel();
p.pack('h',{1/2 []}); % Split it into two columns.
p(1).pack(2,1); % Split the left column into two.
p(2).pack(2,1); % Split the right column into two.
%p.select('all'); % Declare all axes.
%p.identify(); % Show their indices.
ax = gobjects(4,1);
nn = 0;
for ii = 1:2
    for jj = 1:2
        nn = nn+1;
        ax(nn) = p(jj,ii,1).select();
        hold on;
        box on;
    end
end
L = gobjects(4,1);

howmanykb = 10; % The binsize in kb for PosLag figure. If 1, then it's the usual 1kb bins, if 200, then 200 kb bins.
binsize = howmanykb*1000;

mediums = {'LB','MM'};

for mm = 1:2
    medium = mediums{mm};

    if isequal(medium,'LB')
        taudel = 4;
        tau1min = 0;
        tau1max = 30;
    elseif isequal(medium,'MM')
        taudel = 4;
        tau1min = 0;
        tau1max = 44;
    end


    data0 = load('data_shifted.mat').data;
    kGmins = data0.(medium).kGmins;
    kG = data0.(medium).kG;

    % Temporarily load in variables.
    x1 = data0.(medium).oriWT.x';
    logy1 = log(data0.(medium).oriWT.y');

    % To make the fit be a single slope for all data, including both chrs.
    xfull = x1;
    yfull = logy1;

    xpts = [min(x1),0,max(x1)];

    [param,param_err,logyfit] = Fit_SingleSlope_2chr(xfull,yfull,xpts,[],debug);
    logy1fit = logyfit(1:numel(x1));

    forkvel = -kG/param(3);
    forkvel_err = param_err(3)/param(3)*abs(forkvel);
    disp(['Fitted fork velocity: (',num2str(abs(forkvel),3),' ',char(177),' ',num2str(forkvel_err,1),') bp/s']);

    % Normalize with respect to fit.
    if do_ymaxnorm == 1
        logymax = max(logy1fit);
        data.logyfit.full = logy1fit-logymax;
        data.logy.full = logy1-logymax;
    else
        logymax = max(logy1fit);
        logymin = min(logy1fit);
        data.logyfit.full = logy1fit-logymin;
        data.logy.full = logy1-logymin;
    end

    % Get xdel and crtStau from taudel and logyfit.
    xdel = taudel*(x1(2)-x1(1))/(kGmins^(-1)*(logy1fit(2)-logy1fit(1)));

    clear x1 x2 logy1 logy2 xfull yfull xpts1 xpts2 xpts logyfit logy1fit logy2fit forkvel forkvel_err param param_err;

    LRs = {'L','R'};
    fields = {'x','y','logy','logyfit','tau'};

    % Load data into a nicer data structure for plotting.
    data.x.full = data0.(medium).oriWT.x';
    data.y.full = exp(data.logy.full);
    if do_ymaxnorm == 1
        data.tau.full = -kGmins^(-1)*data.logy.full;
    else
        data.tau.full = -kGmins^(-1)*(data.logy.full+logymin-logymax);
    end
    data.sig = data0.(medium).oriWT.sigma_logy;

    for lr = 1:2
        LR = LRs{lr};
        if lr == 1
            xind = (data.x.full < 0);
        else
            xind = (data.x.full >= 0);
        end

        for ii = 1:numel(fields)
            data.(fields{ii}).(LR) = data.(fields{ii}).full(xind);
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if isequal(medium,'LB')
        ii = 0;
    else
        ii = 2;
    end
    sigma = data.sig;
    for lr = 1:2
        LR = LRs{lr};
        ii = ii + 1;

        x = data.x.(LR);
        y = data.y.(LR);
        logy = data.logy.(LR);
        logyfit = data.logyfit.(LR);
        tau = data.tau.(LR);

        taupts = tau1min:taudel:tau1max;

        taumids = mean([taupts(1:end-1);taupts(2:end)],1);

        xpts = sort((-1)^(lr)*(0:xdel:max(abs(x))));
        xmids = mean([xpts(1:end-1);xpts(2:end)],1);

        ypts = ones(size(xpts));
        guess = [ypts,xpts];

        [ypts,ypts_err,~,rep,rep_err] = Fit_PiecewiseLinear_nLeg_SetX(x,logy,guess,kG,sigma,debug);
        rep = abs(rep);
        rep_err = abs(rep_err);

        tmp_rep = nan(size(y));
        tmp_rep(1:ceil(numel(tmp_rep)/numel(rep)):end) = rep;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Top left: logN vs. L. Full L.

        axes(ax(1));

        x_tmp = x;
        xpts_tmp = xpts;

        scatter(x_tmp,y,1,colors{ii},'filled','MarkerFaceAlpha',1);

        if do_ypts == 1
            if do_yptserrbar == 1
                errorbar(xpts_tmp,exp(ypts),exp(ypts).*ypts_err./ypts,'k.-','CapSize',3);
            else
                plot(xpts_tmp,exp(ypts),'k.-');
            end
        end
        %plot(x_tmp,exp(logyfit),'Color','k','LineWidth',1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Bottom left: v vs. L, sharing x axis with top left.

        axes(ax(3));
        L(ii) = plot(xmids,rep,marker,'Color',colors{ii},'LineWidth',0.5,'MarkerSize',12);
        patch([xmids,fliplr(xmids)], [rep-rep_err,fliplr(rep+rep_err)], colors{ii}, 'FaceAlpha',0.3, 'EdgeColor','none');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Top right: Distance vs. time. Lag time starts at 0.

        axes(ax(2));
        if do_varbinsize == 1
            x_var_edges = floor(min(x)/binsize)*binsize:binsize:ceil(max(x)/binsize)*binsize;
            x_var = mean([x_var_edges(1:end-1);x_var_edges(2:end)],1);
            y_var = zeros(size(x_var));
            y_var_sigma = zeros(size(x_var));
            for jj = 1:numel(x_var_edges)-1
                left = x_var_edges(jj);
                right = x_var_edges(jj+1);
                % The error for each of these binned points is just the standard
                % deviation of the individual 1kbp bin points. Otherwise, there's
                % no error and point.
                if sum(x >= left & x < right) == 1
                    y_var(jj) = NaN;
                    y_var_sigma(jj) = NaN;
                else
                    y_var(jj) = mean(y(x >= left & x < right),'omitnan');
                    y_var_sigma(jj) = std(y(x >= left & x < right),'omitnan');
                end
            end
            nnanind_var = ~isnan(y_var);
            x_var = (-1)^lr*x_var(nnanind_var);
            y_var = y_var(nnanind_var);
            sig_var = y_var_sigma(nnanind_var);
            logy_var = log(y_var);
            tau_var = -kGmins^(-1)*logy_var;
            tausig_var = -kGmins^(-1)*sig_var;

            patch([tau_var-tausig_var,fliplr(tau_var+tausig_var)],[x_var,fliplr(x_var)],...
                colors{ii},'FaceAlpha',0.5,'EdgeColor','none');
        else
            if lr == 1
                tauL = tau;
                xL = -x;
            else
                tauR = tau;
                xR = x;
                sizediff = numel(tauL)-numel(tauR);
                if sizediff >= 0
                    taus = [tauL;tauR,nan(1,sizediff)];
                    xs = [xL;xR,nan(1,sizediff)];
                else
                    taus = [tauL,nan(1,-sizediff);tauR];
                    xs = [xL,nan(1,-sizediff);xR];
                end
                rr = logical(randi([0,1],1,length(taus)));
                msize = 3;
                plot(taus(1,rr),xs(1,rr),'.','Color',colors{ii-1},'MarkerSize',msize);
                plot(taus(2,rr),xs(2,rr),'.','Color',colors{ii},'MarkerSize',msize);
                plot(taus(2,~rr),xs(2,~rr),'.','Color',colors{ii},'MarkerSize',msize);
                plot(taus(1,~rr),xs(1,~rr),'.','Color',colors{ii-1},'MarkerSize',msize);

            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Bottom right: v vs. time, sharing x axis with top right.

        axes(ax(4));
        if isequal(LR,'L')
            taumids = fliplr(taumids);
        end
        plot(taumids,rep,marker,'Color',colors{ii},'LineWidth',0.5,'MarkerSize',12);
        patch([taumids,fliplr(taumids)], [rep-rep_err,fliplr(rep+rep_err)], colors{ii}, 'FaceAlpha',0.3, 'EdgeColor','none');

    end

end

% Formatting:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Top left: logN+lag time vs. L. Full L.

axes(ax(1));

% Removing xticklabel.
ax(1).XTickLabel = [];

% Setting y limits.
if do_ymaxnorm == 1
    YLimL = [0.3,1.1];
else
    YLimL = [0.8,3.2];
end
ylim(YLimL);

% Log scale.
set(gca,'YScale','log')

% Include lag time axis on the right.
%{
yyaxis right
ylabel('Lag time: \tau_L (min)');
set(get(gca,'YLabel'),'Rotation',-90,'VerticalAlignment','Bottom')
YLimR = -kGmins^-1*log(YLimL(end:-1:1));
ylim(YLimR)
set( gca, 'ydir', 'reverse' );
yyaxis left;
%}

% oriC line.
xline(0,':');
text(0,YLimL(2),'{\itoriC}','HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',9);

ylabel(ax(1),'Number: N');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bottom left: v vs. L, sharing x axis with top left.

axes(ax(3));

% Setting limits.
ylim([0,3000]);
set(p(1).de.axis,'XLim',[-2.5e6-0.1e6,2.5e6+0.1e6]);

% Making ytick labels be in kb/s.
yts = 0:500:3000;
ytls = yts/1000;
set(ax(3),'YTick',yts,'YTickLabel',ytls);

% oriC line.
xline(0,':');

% xlabel and ylabel.
xlabel(ax(3),'Position: x (bp)');
ylabel(ax(3),'Fork velocity: v (kb/s)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Top right: Distance vs. time. Lag time starts at 0.

axes(ax(2));

% Setting ticks.
ylim([-0.5e6,3e6]);
yts = 0:0.5e6:2.5e6;
ytls = yts/1e6;
set(ax(2),'YTick',yts,'YTickLabel',ytls,'XTickLabel',[]);
%ax(2).YAxis.Exponent = 6;

% ori line
yline(0,':');
xline(0,':')

ylabel(ax(2),'Distance: $|\ell|$ (Mb)','Interpreter','latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bottom right: v vs. time, sharing x axis with top right.

axes(ax(4));

% Making ytick labels be in kb/s.
ylim([0,3000]);
yts = 0:500:3000;
ytls = yts/1000;
set(ax(4),'YTick',yts,'YTickLabel',ytls);

set(p(2).de.axis,'XLim',[0-3,tau1max+5]);

xlabel('Lag time: \tau (min)');
ylabel(ax(4),'Fork velocity: v (kb/s)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overall figure spacing.

p.de.margin = 2;
p(1).marginright = 25;
p.marginright = 10;

p.fontsize = 10;

axes(ax(2));

if do_legend == 1
    p.marginbottom = 17;
    axes(ax(3));
    legendnames = {'LB L','LB R'};
    lg = legend(L,legendnames);
    lg.Orientation = 'horizontal';
    lg.FontSize = 8;
    lgpos = lg.Position; % Left, bottom, width, height.
    lg.Position(1:2) = [0.5-lgpos(3)/2,0.01]; % Left and bottom
    legend('boxon');
end

axpos = ax(3).Position;
if do_arrows == 1
    xlimtot = ax(3).XLim;
    for LR = 1:2
        oridist = 0.1e6;
        arr_L = 0.5e6;
        ha = annotation('arrow');
        ha.HeadWidth = 8;
        ha.Parent = ax(3);
        ha.X = (-1)^(LR)*[oridist,oridist+arr_L];
        ha.Y = [1800,1800];
        hl = annotation('line');
        hl.Parent = ax(3);
        hl.X = (-1)^(LR)*[oridist,oridist];
        hl.Y = [1750,1850];
    end
end

ABCD = 'ACBD';
if do_ABCD == 1
    for ii = 1:4
        axes(ax(ii));
        axpos = get(gca,'Position');
        tt = text(-0.2,1,['{\bf',ABCD(ii),'}'],'HorizontalAlignment','left',...
            'VerticalAlignment','middle','FontSize',10,'Units','normalized');
        %tt.Position(1:2) = [axpos(1),axpos(2)+axpos(4)];
    end
end

p.export([filename,'.pdf']);
%}


%%








