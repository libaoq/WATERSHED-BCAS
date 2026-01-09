function savaProcessSimpleSelectedPoints(pts_size)

%%
flagDone = 0; % when you are happy with selected/removed points, you get out of loop
flagFirstLoop = 1;
while ~flagDone,
    % cretae slash/backslash depending on PC/UNIX 
    if ispc,
        bsl = '\';
    else
        bsl = '/';
    end;

    currDir = pwd;
    % get the file name of the binned data txt file
    pause(0);
    [binFname, binDirPath] = uigetfile('binned_data.mat', 'Select the file name of the binned data bin file...');
    pause(0);
    cd(binDirPath);
    binned_data = load([binDirPath binFname],'-mat');
    binned_data = binned_data.binned_data_combined;

    % load image parameters...
    image_parameters = load('image_parameters.txt');
    x_pixels = image_parameters(1);  % # X pixels in frame
    dataStruct.x_pixels = x_pixels;
    y_pixels = image_parameters(2);  % # Y pixels
    dataStruct.y_pixels = y_pixels;
    excitationDuration_us = image_parameters(3); % duration of EOM opening (excitation) in us
    dataStruct.excitationDuration_us = excitationDuration_us;
    collectionDuration_us = image_parameters(4);     % duration of total collection time (EOM open + EOM closed) in us
    dataStruct.collectionDuration_us = collectionDuration_us;
    numSinglePointAverages = image_parameters(5);            % how many times each pixel was averaged
    dataStruct.numSinglePointAverages = numSinglePointAverages;
    %numSinglePointAverages = uint32(numSinglePointAverages);

    t = binned_data (1,:);
    dataStruct.collectionDurationTimePoints_us = t;
    binned_data = binned_data (2:end,:);
    [numTotReps numTimePtsCollectionDuration] = size(binned_data);
    dataStruct.binned_data = binned_data;
    dataStruct.numTotReps = numTotReps;
    total_binned_pixels = numTotReps;

    bin_interval_us = t(2);
    dataStruct.binInterval_us = bin_interval_us;

    %decay_times = ((excitationDuration_us/bin_interval_us):1:(collectionDuration_us/bin_interval_us));
    %decay_times = ceil(decay_times);
    %dataStruct.decayTimes = decay_times;

    selected_points_indices = load('selected_points.txt');
    dataStruct.selectedPointsEOMtimeIndices = selected_points_indices;
    numImagedPoints = size(selected_points_indices,1);
    dataStruct.numImagedPoints = numImagedPoints;    

% Sergei's 2P-pO2 dye fit parameters pO2 = y0+A1exp(-tau/t1)+A2exp(-tau/t2)
    %sy0 = -9.37638;
    %sA1 = 5686.40211;
    %st1 = 3.58341;
    %sA2 = 269.12134;
    %st2 = 14.52748;
    
% Sergei's 2P-pO2 dye fit (April 2010) parameters pO2 = y0+A1exp(-tau/t1)+A2exp(-tau/t2)
    %sy0 = -8.08577;
    %sA1 = 14768.37256;
    %st1 = 2.22079;
    %sA2 = 237.76688;
    %st2 = 12.38812;
    
% Sergei's 2P-pO2 dye fit (August 05 2010) parameters pO2 = y0+A1exp(-tau/t1)+A2exp(-tau/t2)
    %sy0 = -8.22035;
    %sA1 = 21414.97802;
    %st1 = 2.79095;
    %sA2 = 357.38082;
    %st2 = 11.06016;
    
% Sergei's 2P-pO2 dye fit (Sep 02 2010) parameters pO2 = y0+A1exp(-tau/t1)+A2exp(-tau/t2)
    %sy0 = -9.37637;
    %sA1 = 9510.23021;
    %st1 = 3.20179;
    %sA2 = 256.42021;
    %st2 = 13.58042;
    
% Sergei's New Dye: PtG2P (Feb 2016). parameters pO2 = y0+A1exp(-tau/t1)+A2exp(-tau/t2)
    %sy0 = -14.4221867124227;
    %sA1 = 212.750749735243;
    %st1 = 15.1104663727398; 
    %sA2 = 1003.97151464114;
    %st2 = 3.6573;

% Sergei's New Dye: PtG2P (Oct 2016). parameters pO2 = y0+A1exp(-tau/t1)+A2exp(-tau/t2)
% Pulse: 10 us; Delay: 5 us; Collection: 300 us
    sy0 = -18.8533458785921; 
    sA1 = 166.062669454885; 
    st1 = 17.7515658726037; 
    sA2 = 585.449747064441; 
    st2 = 4.21951599885699; 

    % fit function definitions...
    c0 = [1 40 0];
    decay_profile = @(c,xdata) (c(1)*exp(-xdata/c(2))+c(3));
    decay_profile_error = @(c,xdata,ydata) ((ydata - (c(1)*exp(-xdata/c(2))+c(3)))  );

    %MyOptions = optimset('Display','off','LevenbergMarquardt','on');
    MyOptions = optimoptions(@lsqnonlin,'Display','off','Algorithm','levenberg-marquardt');

    upper_coeffs = [1e8 230 100];
    lower_coeffs = [0.0 0.0 0.0];

    %parameters to compute SO2 from pO2

    %Kelman, GR, J Appl Physiol. 21(4): 1375-1376. 1966.
    %Severinghaus, JW J Appl. Physiol.: Respir. Environ. Exercise Physiol. 46(3): 599-602, 1979
    a1 = -8.5322289e3;
    a2 = 2.1214010e3;
    a3 = -6.7073989e1;
    a4 = 9.3596087e5;
    a5 = -3.1346258e4;
    a6 = 2.3961674e3;
    a7 = -6.7104406e1;

    % input delay for processing decays after closing EOM 
    delay_def = {'5'}; % Sergei's setup is dealing with 20 us delay after EOM is closed due to system time resolution...
    if flagFirstLoop,
        delay_input = inputdlg('process decay after delay (us):','input time delay of decay processing ',1,delay_def);
        delayInterval = str2num(delay_input{1});           
    end;
    dataStruct.delayDecayProcessing_us = delayInterval;


    fittedDecayTimePoints_us = 0:bin_interval_us:(collectionDuration_us - excitationDuration_us - mod((collectionDuration_us - excitationDuration_us),bin_interval_us) - bin_interval_us - delayInterval); %make 20 us delay as in Sergei's calibration system
    fittedDecayTimePoints_indices = ceil((excitationDuration_us+delayInterval)/bin_interval_us)+1:1:(collectionDuration_us-bin_interval_us)/bin_interval_us+1;
    dataStruct.fittedDecayTimePoints_us = fittedDecayTimePoints_us;
    dataStruct.fittedDecayTimePoints_indices = fittedDecayTimePoints_indices;  

    % reserve memory for processed data (from selected sub decays only)
    mean_lifetime_data       = zeros(numImagedPoints, numTimePtsCollectionDuration);
    summed_lifetime_data     = zeros(numImagedPoints, numTimePtsCollectionDuration);
    normalized_lifetime_data = zeros(numImagedPoints, numTimePtsCollectionDuration);
    stdev_lifetime_data      = zeros(numImagedPoints, numTimePtsCollectionDuration);

    lifetimes   =  zeros(numImagedPoints,3);
    pO2         =  zeros(numImagedPoints,1);
    pO2_std_err =  zeros(numImagedPoints,1);
    pO2_rel_err =  zeros(numImagedPoints,1);
    lt_confidence_intervals  = zeros(numImagedPoints, 2);
    lt_std_err  =  zeros(numImagedPoints, 1);

    lt_fit_time =  0:(collectionDuration_us - excitationDuration_us - delayInterval -1);
    lt_fit_profile = zeros(numImagedPoints,length(lt_fit_time));

    % reserve memory for all data (from combi9ned selected and excluded sub
    % decays)
    mean_lifetime_data_all       = zeros(numImagedPoints, numTimePtsCollectionDuration);
    summed_lifetime_data_all     = zeros(numImagedPoints, numTimePtsCollectionDuration);
    normalized_lifetime_data_all = zeros(numImagedPoints, numTimePtsCollectionDuration);
    stdev_lifetime_data_all      = zeros(numImagedPoints, numTimePtsCollectionDuration);

    lifetimes_all   =  zeros(numImagedPoints,3);
    pO2_all         =  zeros(numImagedPoints,1);
    pO2_std_err_all =  zeros(numImagedPoints,1);
    lt_confidence_intervals_all  = zeros(numImagedPoints, 2);
    lt_std_err_all  =  zeros(numImagedPoints, 1);

    lt_fit_profile_all = zeros(numImagedPoints,length(lt_fit_time));

    % input how many decays to average for a single lifetime fit
    % for example, if total # decays is 1000 and we average 100 decays for one
    % fit, this produces 10 lifetimes to average into final result...

    numRepsInSubDecay_def = {num2str(numSinglePointAverages)};  %default is total number of decays
    if flagFirstLoop,
        numRepsInSubDecay_input = inputdlg(['brake total num. decays of ' num2str(numSinglePointAverages) 'into groups of: '],'input number sub-decays for one lifetime ',1,numRepsInSubDecay_def);
        numRepsInSubDecay = str2num(numRepsInSubDecay_input{1});           
    end;
    dataStruct.numRepsInSubDecay = numRepsInSubDecay;

    numSubDecayGroups = round(double(numSinglePointAverages)/numRepsInSubDecay);

    if flagFirstLoop,
        excludePoints = zeros(numImagedPoints,1); % 0 to keep the point, 1 to exclude the point
        excludeSubGroups = zeros(numImagedPoints,numSubDecayGroups); % 0 is subDecayGroup is o.k., 1 if needs to be excluded
    end;
    

    dataStruct.numSubDecayGroupsnumSubDecayGroupsnumSubDecayGroups = numSubDecayGroups;

    % reserve memory for sub decay processing
    sub_summed_lifetime_data  = zeros (numImagedPoints, numSubDecayGroups, numTimePtsCollectionDuration);
    sub_mean_lifetime_data  = zeros (numImagedPoints, numSubDecayGroups, numTimePtsCollectionDuration);
    sub_stdev_lifetime_data  = zeros (numImagedPoints, numSubDecayGroups, numTimePtsCollectionDuration);
    sub_normalized_lifetime_data  = zeros (numImagedPoints, numSubDecayGroups, numTimePtsCollectionDuration);

    sub_lifetimes = zeros(numImagedPoints,numSubDecayGroups,3);
    sub_pO2 = zeros(numImagedPoints,numSubDecayGroups);
    sub_pO2_std_err = zeros(numImagedPoints,numSubDecayGroups);
    sub_lt_confidence_intervals  = zeros(numImagedPoints,numSubDecayGroups, 2);
    sub_lt_std_err  = zeros(numImagedPoints, numSubDecayGroups);

    sub_lt_fit_profile = zeros(numImagedPoints,numSubDecayGroups,length(lt_fit_time));

    subDecIdx = 1:numRepsInSubDecay;
    for j = 1:numImagedPoints,

        disp(['Processing point ' num2str(j) ' out of ' num2str(numImagedPoints) ' points...']); pause(0.001);
        decay_index_offset = (j-1)*numSinglePointAverages;
        decay_indices_all = [1:numSinglePointAverages] + decay_index_offset;

        decay_indices = [];
        for kkk = 1:numSubDecayGroups, % create indexes only for sub decays that are acceptable
            if excludeSubGroups(j,kkk)==0,
                decay_indices = [decay_indices (subDecIdx+(kkk-1)*numRepsInSubDecay)];
            end;
        end;  
        decay_indices = decay_indices+decay_index_offset;

        if isempty(decay_indices), % if all sub groups are excluded, this point is not valid
            excludePoints(j)=1;
        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % remove outliers
            bbb = sum(binned_data(decay_indices,:),2);
            [bbbidx,bbboutliers] = deleteOutliers(bbb);
            decay_indices(bbbidx) = [];
            removedOutliers(j) = length(bbbidx);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            summed_lifetime_data (j,:) = sum(binned_data(decay_indices,:),1);
            mean_lifetime_data (j,:) = mean(binned_data(decay_indices,:),1);
            stdev_lifetime_data (j,:) = std(binned_data(decay_indices,:),1);
            if summed_lifetime_data (j,fittedDecayTimePoints_indices(1)) == 0,
                normalized_lifetime_data (j,:) = 0;
            else
                normalized_lifetime_data (j,:) = summed_lifetime_data (j,:)./summed_lifetime_data (j,fittedDecayTimePoints_indices(1));
            end;
            c0(1) = normalized_lifetime_data(j,1);
            c0(2) = squeeze(sum(normalized_lifetime_data(j,fittedDecayTimePoints_indices))) * bin_interval_us / c0(1);
            c0(3) = 0;

            [cAll,resnorm,resid,exitflag,output,lambda,jacobian_c] = ...
                lsqnonlin(decay_profile_error, c0, [],[],MyOptions,fittedDecayTimePoints_us',...
                                   (normalized_lifetime_data(j,fittedDecayTimePoints_indices))');
            %lsqnonlin(decay_profile_error, c0, lower_coeffs,upper_coeffs,MyOptions,fittedDecayTimePoints_us',(normalized_lifetime_data(j,fittedDecayTimePoints_indices))');
            lifetimes(j,:) = cAll;

            c_ci = nlparci(cAll,resid,'jacobian',jacobian_c,'alpha',0.05);

            lt_confidence_intervals (j,:) = c_ci(2,:);
            lt_std_err(j) = (c_ci(2,2) - c_ci(2,1)) / 3.92;   % assume the CI is symmetric, 95% corresponds to 1.96 * std err

            pO2(j) = sA1.*exp(-lifetimes(j,2)./st1)+sA2.*exp(-lifetimes(j,2)./st2)+sy0;
            pO2_std_err(j) = sqrt((     (-(sA1/st1).*exp(-lifetimes(j,2)./st1)*lt_std_err(j)-(sA2/st2).*exp(-lifetimes(j,2)./st2)*lt_std_err(j) )    )^2);
            pO2_rel_err(j) = pO2_std_err(j) / pO2(j);

            lt_fit_profile (j,:) = feval(decay_profile, cAll, lt_fit_time);
            %             legend_names(j) =  cellstr([num2str(j,2) ' : ' num2str(lifetimes(j,2),4) ' \mus']);
            legend_names(j) =  cellstr([num2str(j,2) ' : ' num2str(pO2(j),3) ' \pm ' num2str(pO2_std_err(j),3) ' mmHg']);
        end;

        % do the same processing on data from all sub decays (accepted and
        % excluded)
        summed_lifetime_data_all (j,:) = sum(binned_data(decay_indices_all,:),1);
        mean_lifetime_data_all (j,:) = mean(binned_data(decay_indices_all,:),1);
        stdev_lifetime_data_all (j,:) = std(binned_data(decay_indices_all,:),1);
        if summed_lifetime_data_all (j,fittedDecayTimePoints_indices(1)) == 0,
            normalized_lifetime_data_all (j,:) = 0;
        else
            normalized_lifetime_data_all (j,:) = summed_lifetime_data_all (j,:)./summed_lifetime_data_all (j,fittedDecayTimePoints_indices(1));
        end;
        c0(1) = normalized_lifetime_data_all(j,1);
        c0(2) = squeeze(sum(normalized_lifetime_data_all(j,fittedDecayTimePoints_indices))) * bin_interval_us / c0(1);
        c0(3) = 0;

        [cAll,resnorm,resid,exitflag,output,lambda,jacobian_c] = ...
            lsqnonlin(decay_profile_error, c0, [],[],MyOptions,fittedDecayTimePoints_us',...
                               (normalized_lifetime_data_all(j,fittedDecayTimePoints_indices))');
        %lsqnonlin(decay_profile_error, c0, lower_coeffs,upper_coeffs,MyOptions,fittedDecayTimePoints_us',(normalized_lifetime_data(j,fittedDecayTimePoints_indices))');
        lifetimes_all(j,:) = cAll;

        c_ci = nlparci(cAll,resid,'jacobian',jacobian_c,'alpha',0.05);

        lt_confidence_intervals_all (j,:) = c_ci(2,:);
        lt_std_err_all(j) = (c_ci(2,2) - c_ci(2,1)) / 3.92;   % assume the CI is symmetric, 95% corresponds to 1.96 * std err

        pO2_all(j) = sA1.*exp(-lifetimes_all(j,2)./st1)+sA2.*exp(-lifetimes_all(j,2)./st2)+sy0;
        pO2_std_err_all(j) = sqrt((     (-(sA1/st1).*exp(-lifetimes_all(j,2)./st1)*lt_std_err_all(j)-(sA2/st2).*exp(-lifetimes_all(j,2)./st2)*lt_std_err_all(j) )  *lt_std_err_all(j))^2);
        
        
        lt_fit_profile_all (j,:) = feval(decay_profile, cAll, lt_fit_time);
        %             legend_names(j) =  cellstr([num2str(j,2) ' : ' num2str(lifetimes(j,2),4) ' \mus']);
        legend_names(j) =  cellstr([num2str(j,2) ' : ' num2str(pO2_all(j),3) ' \pm ' num2str(pO2_std_err_all(j),3) ' mmHg']);



        % now do calculation in subdecays...
        for idxSubDecays = 1:numSubDecayGroups,

            subdecay_index_offset = (j-1)*numSinglePointAverages+(idxSubDecays-1)*numRepsInSubDecay;
            sub_decay_indices = [1:numRepsInSubDecay] + subdecay_index_offset;

            sub_summed_lifetime_data (j,idxSubDecays,:) = sum(binned_data(sub_decay_indices,:),1);
            sub_mean_lifetime_data (j,idxSubDecays,:) = mean(binned_data(sub_decay_indices,:),1);
            sub_stdev_lifetime_data (j,idxSubDecays,:) = std(binned_data(sub_decay_indices,:),1);
            if sub_summed_lifetime_data (j,idxSubDecays,fittedDecayTimePoints_indices(1)) == 0,
                sub_normalized_lifetime_data (j,idxSubDecays,:) = 0;
            else
                sub_normalized_lifetime_data (j,idxSubDecays,:) = sub_summed_lifetime_data (j,idxSubDecays,:)./sub_summed_lifetime_data (j,idxSubDecays,fittedDecayTimePoints_indices(1));
            end;

            %c0(1) = sub_normalized_lifetime_data(j,decay_indices(1));
            %c0(2) = sum(squeeze( sub_normalized_lifetime_data(j,fittedDecayTimePoints_indices(1)) )) * bin_interval_us / c0(1);
            %c0(3) = 0;
            c0 = cAll;

            decayIntensities=squeeze(sub_normalized_lifetime_data(j,idxSubDecays,fittedDecayTimePoints_indices));
            [c,resnorm,resid,exitflag,output,lambda,jacobian_c] = ...
                lsqnonlin(decay_profile_error, c0, [],[],MyOptions,fittedDecayTimePoints_us',decayIntensities);
                %lsqnonlin(decay_profile_error, c0, lower_coeffs,upper_coeffs,MyOptions,fittedDecayTimePoints_us',decayIntensities);
            sub_lifetimes(j,idxSubDecays,:) = c;

            c_ci = nlparci(c,resid,'jacobian',jacobian_c,'alpha',0.05);

            sub_lt_confidence_intervals (j,idxSubDecays,:) = c_ci(2,:);
            sub_lt_std_err(j,idxSubDecays) = (c_ci(2,2) - c_ci(2,1)) / 3.92;   %%assume the CI is symmetric, 95% corresponds to 1.96 * std err


            sub_pO2(j,idxSubDecays) = sA1.*exp(-sub_lifetimes(j,idxSubDecays,2)./st1)+sA2.*exp(-sub_lifetimes(j,idxSubDecays,2)./st2)+sy0;
            sub_pO2_std_err(j,idxSubDecays) = sqrt((     (-(sA1/st1).*exp(-sub_lifetimes(j,idxSubDecays,2)./st1)*sub_lt_std_err(j,idxSubDecays)...
                -(sA2/st2).*exp(-sub_lifetimes(j,idxSubDecays,2)./st2)*sub_lt_std_err(j,idxSubDecays) )  *sub_lt_std_err(j,idxSubDecays))^2);

            sub_lt_fit_profile (j,idxSubDecays,:) = feval(decay_profile, c, lt_fit_time);
            %             legend_names(j) =  cellstr([num2str(j,2) ' : ' num2str(lifetimes(j,2),4) ' \mus']);
            sub_legend_names(j,idxSubDecays) =  cellstr([num2str(j,2) ', ' num2str(idxSubDecays,3) ' : ' num2str(sub_pO2(j,idxSubDecays),3) ' \pm ' num2str(sub_pO2_std_err(j,idxSubDecays),3) ' mmHg']);

        end;
    end;

    stdev_lifetime_data = stdev_lifetime_data./ sqrt(double(numSinglePointAverages));

    dataStruct.sub_mean_lifetime_data = sub_mean_lifetime_data;
    dataStruct.sub_summed_lifetime_data = sub_summed_lifetime_data;
    dataStruct.sub_normalized_lifetime_data = sub_normalized_lifetime_data;
    dataStruct.sub_stdev_lifetime_data = sub_stdev_lifetime_data;
    dataStruct.sub_lifetimes = sub_lifetimes;
    dataStruct.sub_pO2 = sub_pO2;
    dataStruct.sub_pO2_std_err = sub_pO2_std_err;
    dataStruct.sub_lt_confidence_intervals = sub_lt_confidence_intervals;
    dataStruct.sub_lt_std_err = sub_lt_std_err;
    dataStruct.sub_lt_fit_profile = sub_lt_fit_profile;

    dataStruct.mean_lifetime_data = mean_lifetime_data;
    dataStruct.summed_lifetime_data = summed_lifetime_data;
    dataStruct.normalized_lifetime_data = normalized_lifetime_data;
    dataStruct.stdev_lifetime_data = stdev_lifetime_data;
    dataStruct.lifetimes = lifetimes;
    dataStruct.pO2 = pO2;
    dataStruct.pO2_std_err = pO2_std_err;
    dataStruct.pO2_rel_err = pO2_rel_err;
    dataStruct.lt_confidence_intervals = lt_confidence_intervals;
    dataStruct.lt_std_err = lt_std_err;
    dataStruct.lt_fit_time = lt_fit_time;
    dataStruct.lt_fit_profile = lt_fit_profile;

    mean_pO2 = zeros(numImagedPoints,1);
    for j=1:numImagedPoints,
        if excludePoints(j)==0,
            mean_pO2(j) = squeeze(mean(  squeeze(sub_pO2(j,find(squeeze(excludeSubGroups(j,:))==0)) )));
        end;
    end;
    dataStruct.mean_pO2 = mean_pO2;

    dataStruct.ALL.mean_lifetime_data_all = mean_lifetime_data_all;
    dataStruct.ALL.summed_lifetime_data_all = summed_lifetime_data_all;
    dataStruct.ALL.normalized_lifetime_data_all = normalized_lifetime_data_all;
    dataStruct.ALL.stdev_lifetime_data_all = stdev_lifetime_data_all;
    dataStruct.ALL.lifetimes_all = lifetimes_all;
    dataStruct.ALL.pO2_all = pO2_all;
    dataStruct.ALL.pO2_std_err_all = pO2_std_err_all;
    dataStruct.ALL.lt_confidence_intervals_all = lt_confidence_intervals_all;
    dataStruct.ALL.lt_std_err_all = lt_std_err_all;
    dataStruct.ALL.lt_fit_time_all = lt_fit_time;
    dataStruct.ALL.lt_fit_profile_all = lt_fit_profile_all;

    mean_pO2_all = zeros(numImagedPoints,1);
    mean_pO2_all = squeeze(mean(sub_pO2,2));
    dataStruct.ALL.mean_pO2_all = mean_pO2_all;



    %         legend_names = num2str((1:numImagedPoints)');



    %%
    % import survey scan image...
    if flagFirstLoop,
        [surveyFilename, survayPathname] = uigetfile({'*.txt';'*.mat'},'Select a lifetime image (.mat) or survey scan raw data image (.txt): '); % choose lifetime or survey scan image
    
        % act depending on which data was loaded - lifetime image or survey
        % scan raw data
        [Xpathstr, Xname, Xext] = fileparts(surveyFilename);
        if strcmp(Xext,'.mat'),
            load ([survayPathname surveyFilename]);
        else
            tmpImg = load([survayPathname surveyFilename],'-ascii');
            survey_image = flipud(rot90(reshape(tmpImg,image_parameters(1),image_parameters(2))));
        end;
    end;
    surveyFig = figure(1);
    clf;
    fgPts = subplot(2,2,1);

    survey_image1 = medfilt2(survey_image, [2 2],'symmetric');
    % plot first image and decide on contrast 
    if flagFirstLoop,
        MinMinS = min(survey_image1(:));
        MaxMaxS = max(survey_image1(:));
    end;
    condition = 1;
    while condition,
        hSrv = imagesc(survey_image1,[MinMinS MaxMaxS]); 
        colormap(gray(64));
        cbSrv = colorbar;
        ncmp = 64;
        
        maxSurv = MaxMaxS; %max(survey_image1(:));
        minSurv = MinMinS; %min(survey_image1(:));
        c_Surv = min(ncmp,  round(   (ncmp-1).*(survey_image1-minSurv)./(maxSurv-minSurv)  )+1  );
        c_Surv = c_Surv;
        set(hSrv,'CDataMapping','Direct');
        set(hSrv,'CData',c_Surv);
        pause(0);
        subplot(fgPts);
        caxis([min(c_Surv(:)) max(c_Surv(:))]);
        pause(0);
        set(cbSrv,'YLim',[1 ncmp]);
        
        newTicks = 1:10:ncmp;
        newTickLabels = newTicks.*(maxSurv-minSurv)./(ncmp-1)+(minSurv*ncmp-maxSurv)/(ncmp-1);
        newTickLabelsStr=num2str(newTickLabels(1),2);
        for idxLab = 2:length(newTickLabels),
            newTickLabelsStr = char(newTickLabelsStr,num2str(newTickLabels(idxLab),2));
        end;
        pause(0);
        set(cbSrv,'YTick',newTicks);
        set(cbSrv,'YTickLabel',newTickLabelsStr);
        
        xlabel('X (pixels)');
        ylabel('Y (pixels)');
        set(get(gca,'YLabel'),'FontSize',12);
        set(get(gca,'YLabel'),'FontWeight','Bold');
        set(get(gca,'XLabel'),'FontSize',12);
        set(get(gca,'XLabel'),'FontWeight','Bold');
        set(gca,'FontSize',12);
        set(gca,'FontWeight','Bold');
        title('Points on survey image');
        set(surveyFig,'Position',[150 150 1200 900]);
        set(gca,'Position',[0.07 0.55 0.3 0.4]); 
        set(gca,'nextplot','replacechildren');

        button3 = questdlg('Are you happy with the image contrast?','Repeat selection or not...','No');
        NotHappyContrastCondition = strcmp(button3,'No');
        if NotHappyContrastCondition,
            answer= inputdlg({'Enter minimum intensity level of survey scan image'},'Enter minimum intensity level:',1,{num2str(MinMinS)});  %    {'0'});
            answer1 = cell2struct(answer, 'number', 1);
            MinMinS = str2num(answer1.number);
            answer= inputdlg({'Enter maximum intensity level of survey scan image'},'Enter maximum intensity level:',1,{num2str(MaxMaxS)});  %    {'0'});
            answer1 = cell2struct(answer, 'number', 1);
            MaxMaxS = str2num(answer1.number);
        else
            condition = 0;
        end;
    end;
    dataStruct.surveyImage = survey_image;
    dataStruct.surveyImage1 = survey_image1;
    dataStruct.surveyImageGrayBounds = [MinMinS MaxMaxS];
    
%%
    summary.Data = zeros(numImagedPoints,6);
    summary.Data(:,1) = [1:numImagedPoints];
    summary.Data(:,2) = lifetimes(:,2);
    summary.Data(:,3) = lt_std_err(:);
    summary.Data(:,4) = mean_pO2(:);
    summary.Data(:,5) = pO2_std_err(:);
    summary.Data(:,6) = (summed_lifetime_data (:,fittedDecayTimePoints_indices(1)));
    summary.Columns = {'Point number';'Lifetime';'Lifetime error';'pO2';'pO2 error';'Amplitude (photons)'};
    dataStruct.summary = summary;
    dataStruct.excludePoints = excludePoints;
    dataStruct.excludeSubGroups = excludeSubGroups;
    
    fnameDataStruct = ['dataStruct_delay' num2str(delayInterval) '.mat'];
    [dataStructFileOut, dataStructPathOut] = uiputfile( '*.mat','Create output file name for data structure...',fnameDataStruct);
    Name = [dataStructPathOut dataStructFileOut];
    save(Name,'dataStruct','-mat','-v7.3');
 %%   

    fgPts = subplot(2,2,1);
    
    
    
    hSrv = imagesc(survey_image1,[MinMinS MaxMaxS]); 
    daspect([1 1 1]);
    colormap([gray(64);jet(64)]); %colormap(gray(64));
    cbSrv = colorbar;
    ncmp = 64;

    maxSurv = MaxMaxS; %max(survey_image1(:));
    minSurv = MinMinS; %min(survey_image1(:));
    c_Surv = min(ncmp,  round(   (ncmp-1).*(survey_image1-minSurv)./(maxSurv-minSurv)  )+1  );
    c_Surv = c_Surv;
    set(hSrv,'CDataMapping','Direct');
    set(hSrv,'CData',c_Surv);
    pause(0);
    subplot(fgPts);
    caxis([min(c_Surv(:)) max(c_Surv(:))]);
    pause(0);
    set(cbSrv,'YLim',[1 ncmp]);

    newTicks = 1:10:ncmp;
    newTickLabels = newTicks.*(maxSurv-minSurv)./(ncmp-1)+(minSurv*ncmp-maxSurv)/(ncmp-1);
    newTickLabelsStr=num2str(newTickLabels(1),2);
    for idxLab = 2:length(newTickLabels),
        newTickLabelsStr = char(newTickLabelsStr,num2str(newTickLabels(idxLab),2));
    end;
    pause(0);
    set(cbSrv,'YTick',newTicks);
    set(cbSrv,'YTickLabel',newTickLabelsStr);
    
    xlabel('X (pixels)');
    ylabel('Y (pixels)');
    set(get(gca,'YLabel'),'FontSize',12);
    set(get(gca,'YLabel'),'FontWeight','Bold');
    set(get(gca,'XLabel'),'FontSize',12);
    set(get(gca,'XLabel'),'FontWeight','Bold');
    set(gca,'FontSize',12);
    set(gca,'FontWeight','Bold');
    title('Points on survey image');
    set(surveyFig,'Position',[150 150 1200 900]);
    set(gca,'Position',[0.07 0.55 0.3 0.4]);
    set(gca,'nextplot','replacechildren');
    % plot position of points on image  
    for i = 1:numImagedPoints,
        if excludePoints(i)==0,
            clr = 'r';
        else
            clr = 'y';
        end;
        rectangle('Position',[(selected_points_indices(i,3) - 1),(selected_points_indices(i,2) - 1), pts_size, pts_size],...
            'Curvature',[1,1],'FaceColor','r');
        daspect ([1,1,1])
        text((selected_points_indices(i,3) - 5),(selected_points_indices(i,2) - 5),[num2str(i)],...
            'FontSize',10,'FontWeight', 'Bold','Color',clr);
    end;
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    % plot pO2 points on second image
    fgPO2 = subplot(2,2,2);
    title('pO2 (mmHg)');
    
    hSurv = imagesc(survey_image1,[MinMinS MaxMaxS]); 
    %colormap gray;
    cbPO2 = colorbar;
    xlabel('X (pixels)');
    ylabel('Y (pixels)');
    set(get(gca,'YLabel'),'FontSize',12);
    set(get(gca,'YLabel'),'FontWeight','Bold');
    set(get(gca,'XLabel'),'FontSize',12);
    set(get(gca,'XLabel'),'FontWeight','Bold');
    set(gca,'FontSize',12);
    set(gca,'FontWeight','Bold');
    title('pO2 values');
    %set(surveyFig,'Position',[150 150 1200 900]);
    set(gca,'Position',[0.55 0.55 0.3 0.4]); 
    set(gca,'nextplot','replacechildren');
    daspect([1 1 1]);
    
    % change colorma
    colormap([gray(64);jet(64)]);
    ncmp=64;

    maxSurv = MaxMaxS; %max(survey_image1(:));
    minSurv = MinMinS; %min(survey_image1(:));
    c_Surv = min(ncmp,  round(   (ncmp-1).*(survey_image1-minSurv)./(maxSurv-minSurv)  )+1  );
    c_Surv = c_Surv;

    set(hSurv,'CDataMapping','Direct');
    set(hSurv,'CData',c_Surv);
    pause(0);
    subplot(fgPO2);
    caxis([min(c_Surv(:)) max(c_Surv(:))]);
    pause(0);
    set(cbPO2,'YLim',[ncmp+1 2*ncmp]);

    % create proper colormap
    cmap1=jet(64);
    CmapN = size(cmap1,1); % length of colormap[N,3] matrix is N
    cmapPO2idx = ones(length(pO2),1); % vector of colormap indexes for pO2
    if flagFirstLoop,
        minPO2 = min(pO2);
        maxPO2 = max(pO2);
    end;
    
    condition = 1;
    while condition
    
        pO2forPlot = pO2;
        pO2forPlot(find(pO2<=minPO2))=minPO2;
        pO2forPlot(find(pO2>=maxPO2))=maxPO2;

        if maxPO2 == minPO2,
            cmapPO2idx = cmapPO2idx.*round(CmapN/2);
        else
            cmapPO2idx = round(   (  CmapN.*(pO2forPlot-minPO2)+maxPO2-pO2forPlot  )    ./(maxPO2-minPO2) );
            cmapPO2idx(find(cmapPO2idx <= 0)) = 1;
            cmapPO2idx(find(cmapPO2idx > CmapN)) = CmapN;
        end;

        % plot pO2 of points on image  
        for i = 1:numImagedPoints,
            if excludePoints(i)==0,
                clr = 'r';
            else
                clr = 'y';
            end;
            colMap = squeeze(cmap1(cmapPO2idx(i),:));
            rectangle('Position',[(selected_points_indices(i,3) - 1),(selected_points_indices(i,2) - 1), pts_size, pts_size],...
                'Curvature',[1,1],'FaceColor',colMap);
            daspect ([1,1,1])
            %text((selected_points_indices(i,3) - 5),(selected_points_indices(i,2) - 5),[num2str(pO2(i),3) ],...
            %    'FontSize',10,'FontWeight', 'Bold','Color',clr);
        end;

        newTicks = 10:10:ncmp;
        newTickLabels = newTicks.*(maxPO2-minPO2)./(ncmp-1)+(minPO2*ncmp-maxPO2)/(ncmp-1);
        newTickLabelsStr=num2str(newTickLabels(1),2);
        for idxLab = 2:length(newTickLabels),
            newTickLabelsStr = char(newTickLabelsStr,num2str(newTickLabels(idxLab),2));
        end;
        pause(0);
        set(cbPO2,'YTick',newTicks+ncmp);
        set(cbPO2,'YTickLabel',newTickLabelsStr);
        
        button3 = questdlg('Are you happy with the pO2 range?','Repeat selection or not...','No');
        NotHappyPO2range = strcmp(button3,'No');
        if NotHappyPO2range,
            answer= inputdlg({'Enter minimum pO2 level in image'},'Enter minimum pO2 level:',1,{num2str(minPO2)});  %    {'0'});
            answer1 = cell2struct(answer, 'number', 1);
            minPO2 = str2num(answer1.number);
            answer= inputdlg({'Enter maximum pO2 level in image'},'Enter maximum pO2 level:',1,{num2str(maxPO2)});  %    {'0'});
            answer1 = cell2struct(answer, 'number', 1);
            maxPO2 = str2num(answer1.number);
        else
            condition = 0;
        end;
    end;
    
%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
    % then plot rel error of pO2
    % plot pO2 points on second image
    fgPO2relErr = subplot(2,2,3);
    title('pO2_relErr');
    
    hSurv = imagesc(survey_image1,[MinMinS MaxMaxS]); 
    %colormap gray;
    cbPO2err = colorbar;
    xlabel('X (pixels)');
    ylabel('Y (pixels)');
    set(get(gca,'YLabel'),'FontSize',12);
    set(get(gca,'YLabel'),'FontWeight','Bold');
    set(get(gca,'XLabel'),'FontSize',12);
    set(get(gca,'XLabel'),'FontWeight','Bold');
    set(gca,'FontSize',12);
    set(gca,'FontWeight','Bold');
    title('pO2 Relative Error');
    %set(surveyFig,'Position',[150 150 1200 900]);
    set(gca,'Position',[0.07 0.05 0.3 0.4]); 
    set(gca,'nextplot','replacechildren');
    daspect([1 1 1]);
    
    % change colorma
    colormap([gray(64);jet(64)]);
    ncmp=64;

    maxSurv = MaxMaxS; %max(survey_image1(:));
    minSurv = MinMinS; %min(survey_image1(:));
    c_Surv = min(ncmp,  round(   (ncmp-1).*(survey_image1-minSurv)./(maxSurv-minSurv)  )+1  );
    c_Surv = c_Surv;

    set(hSurv,'CDataMapping','Direct');
    set(hSurv,'CData',c_Surv);
    pause(0);
    subplot(fgPO2relErr);
    caxis([min(c_Surv(:)) max(c_Surv(:))]);
    pause(0);
    set(cbPO2err,'YLim',[ncmp+1 2*ncmp]);

    % create proper colormap
    pO2_rel_err = abs(pO2_rel_err);
    cmap1=jet(64);
    CmapN = size(cmap1,1); % length of colormap[N,3] matrix is N
    cmapPO2relErrIdx = ones(length(pO2_rel_err),1); % vector of colormap indexes for pO2
    if flagFirstLoop,
        minPO2relErr = min(pO2_rel_err);
        maxPO2relErr = max(pO2_rel_err);
    end;
    
    condition = 1;
    while condition
    
        pO2_rel_err_forPlot = pO2_rel_err;
        pO2_rel_err_forPlot(find(pO2_rel_err<=minPO2relErr))=minPO2relErr;
        pO2_rel_err_forPlot(find(pO2_rel_err>=maxPO2relErr))=maxPO2relErr;

        if maxPO2relErr == minPO2relErr,
            cmapPO2relErrIdx = cmapPO2relErrIdx.*round(CmapN/2);
        else
            cmapPO2relErrIdx = round(   (  CmapN.*(pO2_rel_err_forPlot-minPO2relErr)+maxPO2relErr-pO2_rel_err_forPlot  )    ./(maxPO2relErr-minPO2relErr) );
            cmapPO2relErrIdx(find(cmapPO2relErrIdx <= 0)) = 1;
            cmapPO2relErrIdx(find(cmapPO2relErrIdx > CmapN)) = CmapN;
        end;

        % plot pO2 of points on image  
        for i = 1:numImagedPoints,
            if excludePoints(i)==0,
                clr = 'r';
            else
                clr = 'y';
            end;
            colMap = squeeze(cmap1(cmapPO2relErrIdx(i),:));
            rectangle('Position',[(selected_points_indices(i,3) - 1),(selected_points_indices(i,2) - 1), pts_size, pts_size],...
                'Curvature',[1,1],'FaceColor',colMap);
            daspect ([1,1,1])
            %text((selected_points_indices(i,3) - 5),(selected_points_indices(i,2) - 5),[num2str(pO2(i),3) ],...
            %    'FontSize',10,'FontWeight', 'Bold','Color',clr);
        end;

        newTicks = 10:10:ncmp;
        newTickLabels = newTicks.*(maxPO2relErr-minPO2relErr)./(ncmp-1)+(minPO2relErr*ncmp-maxPO2relErr)/(ncmp-1);
        newTickLabelsStr=num2str(newTickLabels(1),2);
        for idxLab = 2:length(newTickLabels),
            newTickLabelsStr = char(newTickLabelsStr,num2str(newTickLabels(idxLab),2));
        end;
        pause(0);
        set(cbPO2err,'YTick',newTicks+ncmp);
        set(cbPO2err,'YTickLabel',newTickLabelsStr);
        
        button3 = questdlg('Are you happy with the pO2 rel error range?','Repeat selection or not...','No');
        NotHappyPO2relErrRange = strcmp(button3,'No');
        if NotHappyPO2relErrRange,
            answer= inputdlg({'Enter minimum pO2 rel err level in image'},'Enter minimum pO2 rel err level:',1,{num2str(minPO2relErr)});  %    {'0'});
            answer1 = cell2struct(answer, 'number', 1);
            minPO2relErr = str2num(answer1.number);
            answer= inputdlg({'Enter maximum pO2 rel err level in image'},'Enter maximum pO2 rel err level:',1,{num2str(maxPO2relErr)});  %    {'0'});
            answer1 = cell2struct(answer, 'number', 1);
            maxPO2relErr = str2num(answer1.number);
        else
            condition = 0;
        end;
    end;
  %%  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % plot number of collected photons (baseline rremoved) for each point
 
 % summed_lifetime_data     = zeros(numImagedPoints, numTimePtsCollectionDuration);
 numDetectedPhotons = squeeze(sum(summed_lifetime_data (:,fittedDecayTimePoints_indices),2)) - squeeze(mean( summed_lifetime_data(:,fittedDecayTimePoints_indices(end-20:end)) ,2)) .* length(fittedDecayTimePoints_indices);
 
 fgNumPhot = subplot(2,2,4);
    title('Number of Photons in Decay (no baseline)');
    
    hSurv = imagesc(survey_image1,[MinMinS MaxMaxS]); 
    %colormap gray;
    cbNumPhot = colorbar;
    xlabel('X (pixels)');
    ylabel('Y (pixels)');
    set(get(gca,'YLabel'),'FontSize',12);
    set(get(gca,'YLabel'),'FontWeight','Bold');
    set(get(gca,'XLabel'),'FontSize',12);
    set(get(gca,'XLabel'),'FontWeight','Bold');
    set(gca,'FontSize',12);
    set(gca,'FontWeight','Bold');
    title('Number of Photons in Decay (no baseline)');
    %set(surveyFig,'Position',[150 150 1200 900]);
    set(gca,'Position',[0.55 0.05 0.3 0.4]); 
    set(gca,'nextplot','replacechildren');
    daspect([1 1 1]);
    
    % change colorma
    colormap([gray(64);jet(64)]);
    ncmp=64;

    maxSurv = MaxMaxS; %max(survey_image1(:));
    minSurv = MinMinS; %min(survey_image1(:));
    c_Surv = min(ncmp,  round(   (ncmp-1).*(survey_image1-minSurv)./(maxSurv-minSurv)  )+1  );
    c_Surv = c_Surv;

    set(hSurv,'CDataMapping','Direct');
    set(hSurv,'CData',c_Surv);
    pause(0);
    subplot(fgNumPhot);
    caxis([min(c_Surv(:)) max(c_Surv(:))]);
    pause(0);
    set(cbNumPhot,'YLim',[ncmp+1 2*ncmp]);

    % create proper colormap
    cmap1=jet(64);
    CmapN = size(cmap1,1); % length of colormap[N,3] matrix is N
    cmapNumPhotIdx = ones(length(numDetectedPhotons),1); % vector of colormap indexes for pO2
    if flagFirstLoop,
        minNumPhot = min(numDetectedPhotons);
        maxNumPhot = max(numDetectedPhotons);
    end;
    
    condition = 1;
    while condition
    
        numDetectedPhotons_forPlot = numDetectedPhotons;
        numDetectedPhotons_forPlot(find(numDetectedPhotons<=minNumPhot))=minNumPhot;
        numDetectedPhotons_forPlot(find(numDetectedPhotons>=maxNumPhot))=maxNumPhot;

        if maxNumPhot == minNumPhot,
            cmapNumPhotIdx = cmapNumPhotIdx.*round(CmapN/2);
        else
            cmapNumPhotIdx = round(   (  CmapN.*(numDetectedPhotons_forPlot-minNumPhot)+maxNumPhot-numDetectedPhotons_forPlot  )    ./(maxNumPhot-minNumPhot) );
            cmapNumPhotIdx(find(cmapNumPhotIdx <= 0)) = 1;
            cmapNumPhotIdx(find(cmapNumPhotIdx > CmapN)) = CmapN;
        end;

        % plot pO2 of points on image  
        for i = 1:numImagedPoints,
            if excludePoints(i)==0,
                clr = 'r';
            else
                clr = 'y';
            end;
            colMap = squeeze(cmap1(cmapNumPhotIdx(i),:));
            rectangle('Position',[(selected_points_indices(i,3) - 1),(selected_points_indices(i,2) - 1), pts_size, pts_size],...
                'Curvature',[1,1],'FaceColor',colMap);
            daspect ([1,1,1])
            %text((selected_points_indices(i,3) - 5),(selected_points_indices(i,2) - 5),[num2str(pO2(i),3) ],...
            %    'FontSize',10,'FontWeight', 'Bold','Color',clr);
        end;

        newTicks = 10:10:ncmp;
        newTickLabels = newTicks.*(maxNumPhot-minNumPhot)./(ncmp-1)+(minNumPhot*ncmp-maxNumPhot)/(ncmp-1);
        newTickLabelsStr=num2str(newTickLabels(1),2);
        for idxLab = 2:length(newTickLabels),
            newTickLabelsStr = char(newTickLabelsStr,num2str(newTickLabels(idxLab),2));
        end;
        pause(0);
        set(cbNumPhot,'YTick',newTicks+ncmp);
        set(cbNumPhot,'YTickLabel',newTickLabelsStr);
        
        button3 = questdlg('Are you happy with the nnum photons range?','Repeat selection or not...','No');
        NotHappyNumPhotRange = strcmp(button3,'No');
        if NotHappyNumPhotRange,
            answer= inputdlg({'Enter minimum num phot level in image'},'Enter minimum number photons level:',1,{num2str(minNumPhot)});  %    {'0'});
            answer1 = cell2struct(answer, 'number', 1);
            minNumPhot = str2num(answer1.number);
            answer= inputdlg({'Enter maximum num phot level in image'},'Enter maximum num phot level:',1,{num2str(maxNumPhot)});  %    {'0'});
            answer1 = cell2struct(answer, 'number', 1);
            maxNumPhot = str2num(answer1.number);
        else
            condition = 0;
        end;
    end;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
    fnameJPG = ['pO2_points_delay' num2str(delayInterval) '.jpg'];
    [imgNameOut, imgPathOut] = uiputfile( '*.jpg','Create output file name for survey scan image',fnameJPG);
    print(surveyFig,'-djpeg',imgNameOut,'-r300');


    %%
    % plot lifetimes and pO2's to see how imaging went...
    fnameControl = ['controll_delay' num2str(delayInterval) '.tif'];
    [controllNameOut, controllPathOut] = uiputfile( '*.tif','Create output file name for controll image',fnameControl);
    %%
    flagChange = 0; % is there any change in rejected sub decays
    answer= inputdlg({'How many individuall fits to plot'},'Enter # pooint fits in tiff file:',1,{num2str(numImagedPoints)});  %    {'0'});
    numberFitsToPlot = str2num(answer{1});
    
    close(surveyFig);
    
    for nn = 1:numberFitsToPlot,

        hftmp = figure(100+nn); 
        hold off;
        set(hftmp,'Position',get(0,'ScreenSize'));

%         subplot(2,3,1); % plot lifetimes of sub decays
%         plot(sub_lifetimes(nn,:,2),'.-b'); legendstr = {sprintf('%20s','all lifetimes')};
%         hold on;
%         fooIdx = find(squeeze(excludeSubGroups(nn,:))==1);
%         if ~isempty(fooIdx),
%             plot(fooIdx,sub_lifetimes(nn,fooIdx,2),'xr','MarkerSize',15);
%             legendstr = [legendstr sprintf('%20s','excluded lifetimes')];
%         end;
%         title(['Point ' num2str(nn) ' Lifetimes']);xlabel('sub decay indice');ylabel('Lifetime (us)');
%         legend(legendstr);

%         subplot(2,3,2); % plot po2 values  
%         hold on;
%         plot(ones(numSubDecayGroups,1).*mean_pO2(nn),'g'); legendstr = {sprintf('%20s','mean (select.) pO2s')};
%         plot(ones(numSubDecayGroups,1).*pO2(nn),'r'); legendstr = [legendstr sprintf('%20s','true (sel.) pO2s')]; 
%         plot([1:numSubDecayGroups]',sub_pO2(nn,:),'.-b'); legendstr = [legendstr sprintf('%20s','sub decay pO2s')];
%         fooIdx = find(squeeze(excludeSubGroups(nn,:))==1);
%         if ~isempty(fooIdx),
%             plot(fooIdx,sub_pO2(nn,fooIdx),'xr','MarkerSize',15);
%             %legendstr = [legendstr sprintf('%20s','excluded pO2s')];
%         end;
%         title('pO2');xlabel('sub decay indice');ylabel('pO_2 (mmHg)'); legend(legendstr);

%         subplot(2,3,3); % plot sub decay amplitudes...
%         plot(sub_summed_lifetime_data (nn,:,fittedDecayTimePoints_indices(1)),'.-b');
%         legendstr = {sprintf('%20s','all amplitudes')};
%         hold on;
%         fooIdx = find(squeeze(excludeSubGroups(nn,:))==1);
%         if ~isempty(fooIdx),
%             plot(fooIdx,sub_summed_lifetime_data (nn,fooIdx,fittedDecayTimePoints_indices(1)),'xr','MarkerSize',15);
%             legendstr = [legendstr sprintf('%20s','excluded amplitudes')];
%         end;
%         title('Amplitude');legend(legendstr);

%         subplot(2,3,4); % fitted data intensity integral
%         plot(sum(sub_summed_lifetime_data (nn,:,fittedDecayTimePoints_indices),3),'.-b');
%         hold on; legendstr = {sprintf('%30s','fitted data intensity integral')};
%         fooIdx = find(squeeze(excludeSubGroups(nn,:))==1);
%         if ~isempty(fooIdx),
%             plot(fooIdx,sum(sub_summed_lifetime_data (nn,fooIdx,fittedDecayTimePoints_indices),3),'xr','MarkerSize',15);
%             legendstr = [legendstr sprintf('%30s','excluded intensity integral')];
%         end;
%         title('Intensity integral'); legend(legendstr);

        fooIdx = [];
        hflll = subplot(1,1,1); % plot logaritmic decays
        hold on;
        if excludePoints(nn) == 0, 
            semilogy(t(fittedDecayTimePoints_indices), (lifetimes(nn,1).*exp(-fittedDecayTimePoints_us./lifetimes(nn,2))+lifetimes(nn,3)) .* summed_lifetime_data (nn,fittedDecayTimePoints_indices(1)),'k','LineWidth',2 );
            semilogy(t,summed_lifetime_data(nn,:),'.b');
            ltm = ['lifetime = ' sprintf('%3.1f',lifetimes(nn,2)) ' us'];
            legendstr = {sprintf('%20s',ltm)};
            legendstr = [legendstr sprintf('%20s','sum decays')];
        end;
        semilogy(t,squeeze(sub_summed_lifetime_data(nn,:,:)));
        %legendstr = {sprintf('%20s','phosphor decays')};
        if ~isempty(fooIdx),
            semilogy(t,squeeze(sub_summed_lifetime_data(nn,fooIdx,:)),'k');
            %legendstr = [legendstr sprintf('%20s','excluded decays')];
        end;
        set(hflll,'YScale','log');
        xlabel('Time (us)');
        title(['Sub- and total decay. Removed outliers from total: ' num2str(removedOutliers(nn))]);
        if excludePoints(nn)==0,
            legend(legendstr);
        else
            legend('This point is excluded');
        end;

%         subplot(2,3,6); % fited data integral,baseline removed (# counts per decay)
%         hold on;
%         plot(sum(sub_summed_lifetime_data (nn,:,fittedDecayTimePoints_indices),3) - squeeze(mean(sub_summed_lifetime_data(nn,:,fittedDecayTimePoints_indices(end-20:end)),3)).*length(fittedDecayTimePoints_indices),'.-b' );
%         fooIdx = find(squeeze(excludeSubGroups(nn,:))==1);
%         if ~isempty(fooIdx),
%             plot(fooIdx,sum(sub_summed_lifetime_data (nn,fooIdx,fittedDecayTimePoints_indices),3) - squeeze(mean(sub_summed_lifetime_data(nn,fooIdx,fittedDecayTimePoints_indices(end-20:end)),3)).*length(fittedDecayTimePoints_indices),'xr','MarkerSize',15);
%         end;
%         title('Intensity integral (baseline removed)'); 

        print(hftmp,'temp.tif','-dtiff');
        AAA = imread('temp.tif');
        imwrite(AAA,[controllPathOut controllNameOut],'tiff','WriteMode','append');

        close(hftmp);
        % ask if some sub-decays are to be removed
        fooIdx = find(squeeze(excludeSubGroups(nn,:))==1);
        if isempty(fooIdx),
            excludedSoFar = {'none'};
        else
            excludedSoFar = ' ';
            for ifoo = 1:length(fooIdx),
                excludedSoFar = [excludedSoFar ' ' num2str(fooIdx(ifoo))];
            end;
            excludedSoFar = {excludedSoFar};
        end;
        %excludeInput = inputdlg('Which sub decays to exclude (1 2 5...) or none?','Input sub decays to exclude ',1,excludedSoFar);
        excludeInput{1} = 'none';
        if ~strcmp(excludeInput{1},'none'),
            excludeInput = str2num(excludeInput{1});
            if ~isequal(excludeInput, fooIdx),
                flagChange = 1;
                excludeSubGroups(nn,excludeInput)=1;
                if isempty(find(excludeSubGroups(nn,:)==0)),
                    excludePoints(nn)=1;
                end;
            end;
        else
            if ~isempty(fooIdx),
                flagChange = 1;
                excludeSubGroups(nn,:)=0;
                excludePoints(nn)=0;
            end;
        end;

    end;
    
    if flagChange,
        delete([controllPathOut controllNameOut]);
    end;
    delete('temp.tif');
    close all;
    flagDone = ~flagChange;
    flagFirstLoop = 0; 

end;

