function S = addToTbl_KKN(data,varargin)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % data = ' 

    p = inputParser;
    addRequired(p, 'data',  @(x) (ischar(x)) || isstruct(x));
    addParameter(p, 'SAVE_NAME', [], @ischar);

    parse(p, data, varargin{:});
    data = p.Results.data;
    SAVE_NAME = p.Results.SAVE_NAME;

    if ischar(data) 
        load(data,'S');
        fprintf('\n----Data loaded ----\n')
        if ~isempty(SAVE_NAME)
            [filepath,~,fileext] = fileparts(data);
            table_path = fullfile(filepath, [S.sess_name '-', SAVE_NAME, fileext]);
        end
    else
        S = data;
    end

    fprintf('\n------------------------------\n')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. Edits to clust tables based on sorting
    lastClust = 0; mins = []; maxs = [];
    for c = 1:size(S.kilosort,1)
        clusts = S.kilosort(c).clusters;
        if contains(string(clusts.sess_name(1)),'scrappy')
            clusts.monkey = repmat(categorical("S"),height(clusts),1);
        else
            clusts.monkey = repmat(categorical("Y"),height(clusts),1);
        end

        clusts.row_id = lastClust + clusts.cluster_id + 1;
        lastClust = clusts.row_id(end);

        mins = [mins min(S.kilosort(c).spike_times)];
        maxs = [maxs max(S.kilosort(c).spike_times)];

        S.kilosort(c).clusters = clusts;
    end

    % bw = 10*30000;
    % spks_binned = nan(lastClust,numel(min(mins):bw:max(maxs)-bw));
    % x = (min(mins):bw:max(maxs)-bw)./30000;
    % clust = 1;
    % for c = 1:size(S.kilosort,1)
    %     clusts = S.kilosort(c).clusters;
    %     for u = 1:height(clusts)
    %         spks_binned(clust,:) = bin_spktimes(S.kilosort(c).spike_times(S.kilosort(c).spike_clusters==clusts.cluster_id(u)), 'StartTime', min(mins), 'EndTime', max(maxs), 'BinSize', bw);
    %         clust = clust + 1;
    %     end
    % end
    % mean_fr = mean(spks_binned);
    % idx = findchangepts(mean_fr, 'Statistic', 'std');
    % cutoff_sec = x(idx);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEHAVIOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    fns = fieldnames(S);

    % 2. Make some edits to RFM
    matchingFields = fns(contains(fns, {'rfmp', 'rfMapping'}, 'IgnoreCase', true));
    for f = 1:numel(matchingFields)
        tbl = S.(matchingFields{f}).tbl;

        if contains(string(tbl.sess_name(1)),'scrappy')
            tbl.monkey = repmat(categorical("S"),height(tbl),1);
        else
            tbl.monkey = repmat(categorical("Y"),height(tbl),1);
        end

        tbl.STIM_ON(tbl.result~='CORRECT') = cellfun(@(q) q(1:end-1), tbl.STIM_ON(tbl.result~='CORRECT'), 'uni', 0);
        tbl.STIM_OFF(tbl.result~='CORRECT') = cellfun(@(q) q(1:end-1), tbl.STIM_OFF(tbl.result~='CORRECT'), 'uni', 0);
        tbl.conditions = cellfun(@(q,v) q(1:numel(v)), tbl.conditions, tbl.STIM_ON, 'uni', 0);
        tbl = tbl(cellfun(@(q) ~isempty(q), tbl.conditions, 'uni', 1),:);

        eyePos = cellfun(@(x) filterEyeTraces_EyeLink(x), tbl.eyedata, 'uni', 0);
        [eyeVel, eyeAcc] = cellfun(@(x) calcDerivative_eyeTraces(x), eyePos, 'uni', 0);

        tbl.eyePos = eyePos; tbl.eyeVel = eyeVel; tbl.eyeAcc = eyeAcc;

        saccades = cellfun(@(q) detect_saccades(q), cellfun(@(v,f) v(:,f(1):end), eyeVel, tbl.FIXATE, 'uni', 0), 'uni', 0);
        tbl.saccades = cellfun(@(q,v) num2cell(q+v(1),2), saccades, tbl.FIXATE, 'uni', 0);

        S.(matchingFields{f}).tbl = tbl;
    end

    % Make some edits to MDIR
    matchingFields = fns(contains(fns, {'mdir', 'dirmem'}, 'IgnoreCase', true));
    for f = 1:numel(matchingFields)
        tbl = S.(matchingFields{f}).tbl;

        if contains(string(tbl.sess_name(1)),'scrappy')
            tbl.monkey = repmat(categorical("S"),height(tbl),1);
        else
            tbl.monkey = repmat(categorical("Y"),height(tbl),1);
        end

        tbl = tbl(tbl.result=='CORRECT',:);

        spks1 = cellfun(@(r) cellfun(@(q) numel(q), r, 'uni', 1), tbl.spiketimes_1, 'uni', 0);
        spks1 = cellfun(@(q,v) (q./v)*1000, spks1, num2cell(tbl.END_TRIAL), 'uni', 0);

        if ismember('spiketimes_2',tbl.Properties.VariableNames)
            spks2 = cellfun(@(r) cellfun(@(q) numel(q), r, 'uni', 1), tbl.spiketimes_2, 'uni', 0);
            spks2 = cellfun(@(q,v) (q./v)*1000, spks2, num2cell(tbl.END_TRIAL), 'uni', 0);

            spks = [vertcat(spks1{:}) vertcat(spks2{:})];
        else
            spks = vertcat(spks1{:});
        end

        mask = (spks >= (mean(spks, 1) - 3*std(spks, 0, 1))) & (spks <= (mean(spks, 1) + 3*std(spks, 0, 1)));
        tbl = tbl((sum(mask,2) < (mean(sum(mask,2)) + std(sum(mask,2))*3)) & (sum(mask,2) > (mean(sum(mask,2)) - std(sum(mask,2))*3)), :);

        eyePos = cellfun(@(x) filterEyeTraces_EyeLink(x), tbl.eyedata, 'uni', 0);
        [eyeVel, eyeAcc] = cellfun(@(x) calcDerivative_eyeTraces(x), eyePos, 'uni', 0);

        tbl.eyePos = eyePos; tbl.eyeVel = eyeVel; tbl.eyeAcc = eyeAcc;

        saccades = cellfun(@(q) detect_saccades(q), cellfun(@(v,f) v(:,f(1):end), eyeVel, tbl.FIXATE, 'uni', 0), 'uni', 0);
        tbl.saccades = cellfun(@(q,v) num2cell(q+v(1),2), saccades, tbl.FIXATE, 'uni', 0);

        tbl.saccadeOnset = nan(height(tbl),1);
        for t = 1:height(tbl)
            x = cellfun(@(q) q(1), tbl.saccades{t}, 'uni', 1) - tbl.SACCADE(t);
            x(x>0)=NaN;
            [~,m] = min(abs(x));

            tbl.saccadeOnset(t) = tbl.saccades{t}{m}(1);
            tbl.saccadeOffset(t) = tbl.saccades{t}{m}(2);
        end

        tbl.saccadeLatency = tbl.saccadeOnset - cell2mat(tbl.FIX_OFF);
        
        tbl = tbl(tbl.saccadeLatency>=100,:);

        % Det if first saccade out of fixation window landed in targ win
        inTargets = nan(height(tbl),1); 
        [dThetas,dRhos,dists] = deal(cell(height(tbl),1));
        for t = 1:height(tbl)
            % radial position of eye at s
            [theta_eye, rho_eye] = cart2pol(tbl.eyePos{t}(1,tbl.saccadeOffset(t)+50),tbl.eyePos{t}(2,tbl.saccadeOffset(t)+50));
            rho_targ = tbl.distance(t);
            theta_targ = deg2rad(tbl.angle(t));
            r_window = pix2deg(tbl.params(t).block.targWinRad,tbl.params(t).block.screenDistance,tbl.params(t).block.pixPerCM);

            theta_eye = mod(theta_eye, 2*pi);
            theta_targ = mod(theta_targ, 2*pi);

            % Signed difference: positive = clockwise
            dThetas{t} = rad2deg(- (mod(theta_eye - theta_targ + pi, 2*pi) - pi));
            dRhos{t} = rho_eye-rho_targ;

            % Compute distance using law of cosines
            dist = sqrt(rho_eye.^2 + rho_targ^2 - 2*rho_eye*rho_targ.*cos(theta_eye - theta_targ));
            dists{t} = dist;
            
            % Logical array: true if eye is inside target window
            inTargets(t) = dist <= r_window;  
        end

        tbl.saccadeOffset_dTheta = dThetas;
        tbl.saccadeOffset_dRho = dRhos;
        tbl.saccadeOffset_dist = dists;

        tbl = tbl(inTargets==1,:);

        S.(matchingFields{f}).tbl = tbl;
    end
    
    % Make some edits to PURS
    matchingFields = fns(contains(fns, {'purs', 'pursuit'}, 'IgnoreCase', true));
    for f = 1:numel(matchingFields)
        tbl = S.(matchingFields{f}).tbl;

        if contains(string(tbl.sess_name(1)),'scrappy')
            tbl.monkey = repmat(categorical("S"),height(tbl),1);
        else
            tbl.monkey = repmat(categorical("Y"),height(tbl),1);
        end

        tbl = tbl(tbl.result=='CORRECT' & tbl.jump==-1,:);
    
        spks1 = cellfun(@(r) cellfun(@(q) numel(q), r, 'uni', 1), tbl.spiketimes_1, 'uni', 0);
        spks1 = cellfun(@(q,v) (q./v)*1000, spks1, num2cell(tbl.END_TRIAL), 'uni', 0);
    
        if ismember('spiketimes_2',tbl.Properties.VariableNames)
            spks2 = cellfun(@(r) cellfun(@(q) numel(q), r, 'uni', 1), tbl.spiketimes_2, 'uni', 0);
            spks2 = cellfun(@(q,v) (q./v)*1000, spks2, num2cell(tbl.END_TRIAL), 'uni', 0);
    
            spks = [vertcat(spks1{:}) vertcat(spks2{:})];
        else
            spks = vertcat(spks1{:});
        end
    
        mask = (spks >= (mean(spks, 1) - 3*std(spks, 0, 1))) & (spks <= (mean(spks, 1) + 3*std(spks, 0, 1)));
        tbl = tbl((sum(mask,2) < (mean(sum(mask,2)) + std(sum(mask,2))*3)) & (sum(mask,2) > (mean(sum(mask,2)) - std(sum(mask,2))*3)), :);
    
        eyePos = cellfun(@(x) filterEyeTraces_EyeLink(x), tbl.eyedata, 'uni', 0);
        [eyeVel, eyeAcc] = cellfun(@(x) calcDerivative_eyeTraces(x), eyePos, 'uni', 0);
    
        tbl.eyePos = eyePos; tbl.eyeVel = eyeVel; tbl.eyeAcc = eyeAcc;

        saccades = cellfun(@(q) detect_saccades(q, 'VEL_THRESH', 30, 'ACC_THRESH', 500), cellfun(@(v,f) v(:,f(1):end), eyeVel, num2cell(tbl.FIXATE), 'uni', 0), 'uni', 0);
        tbl.saccades = cellfun(@(q,v) num2cell(q+v(1),2), saccades, num2cell(tbl.FIXATE), 'uni', 0);

        [pursuitOnset, pursuitLatency] = cellfun(@(u,v,w) detect_pursuitOnset(u, v, w, 'PLOT_TRACES', false), tbl.eyeVel, num2cell(tbl.PURSUIT_TARG_ON), num2cell(tbl.pursuitSpeed), 'uni', 1); 
        tbl.pursuitOnset = pursuitOnset;
        tbl.pursuitLatency = pursuitLatency;

        tbl = tbl(tbl.pursuitLatency >= 60 & tbl.pursuitLatency < 150,:);

        csTrials = cellfun(@(u,v) sum(((cellfun(@(q) q(2), u, 'uni', 1) - (v)) >= -25) & ((cellfun(@(q) q(1), u, 'uni', 1) - (v+tbl.params(1).block.crossingTime)) < 150)), tbl.saccades, num2cell(tbl.PURSUIT_TARG_ON), 'uni', 1);

        tbl = tbl(~csTrials,:);

        [~,rVel] = cellfun(@(q) cart2pol(q(1,:),q(2,:)), tbl.eyeVel, 'uni', 0);
        rVel = cellfun(@(u,v) u(v+50:v+200), rVel, num2cell(tbl.pursuitOnset), 'uni', 0);
        tbl = tbl(cellfun(@(q,v) min(q) > 0.3*v, rVel, num2cell(tbl.pursuitSpeed), 'uni', 1),:);

        % f = figure;
        % [~,rVel] = cellfun(@(q) cart2pol(q(1,:),q(2,:)), tbl.eyeVel, 'uni', 0);
        % rVel = cellfun(@(u,v) u(v-100:v+300), rVel, num2cell(tbl.pursuitOnset), 'uni', 0);
        % rVel = cellfun(@(x) filterEyeTraces_EyeLink(x,'CUTOFF_FREQUENCY', 20), rVel, 'uni', 0);
        % 
        % x = [0:400]-100;
        % xline(0,'k--');
        % hold on;
        % for t = 1:height(tbl)
        %     plot(x,rVel{t},'k-')
        % end
        % ylim([0,30]);

        S.(matchingFields{f}).tbl = tbl;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VMI, if dirmem was ran
    fields = fieldnames(S);
    matchingFields1 = fields(contains(fields, {'dirmem', 'mdir'}, 'IgnoreCase', true));

    if ~isempty(matchingFields1)
        fprintf('\n~~CALCULATING MDIR METRICS~~\n');
        Tmdir = []; 
        for mm = 1:numel(matchingFields1)
            Tmdir = [Tmdir; S.(matchingFields1{mm}).tbl];
        end

        for prb = 1:size(S.kilosort,1)
            % Percent trials fired on
            vals = cellfun(@(u) cellfun(@(q) numel(q), u, 'uni', 1)>0, Tmdir.(sprintf('spiketimes_%d',prb)), 'uni', 0);
            S.kilosort(prb).clusters.ratio_mdir_trials = sum(vertcat(vals{:}),1)'./(height(Tmdir));
            
            % delay period FR
            delay_fr = cellfun(@(u,v,w) cellfun(@(q) (sum(q>=u(1) & q<v)./(v-u(1)))*1000, w, 'uni', 0), Tmdir.TARG_ON, num2cell(Tmdir.SACCADE), Tmdir.(sprintf('spiketimes_%d',prb)), 'uni', 0);
            delay_fr = cell2mat(vertcat(delay_fr{:}));

            % dprime between evoked and spotantaneous activity
            vis_fr = cellfun(@(u,w) cellfun(@(q) sum(q>=(u(1)+50) & q<(u(1)+250)), w, 'uni', 0), Tmdir.TARG_ON, Tmdir.(sprintf('spiketimes_%d',prb)), 'uni', 0);
            vis_fr = cell2mat(vertcat(vis_fr{:}));

            sac_fr = cellfun(@(u,w) cellfun(@(q) sum(q>=(u-75) & q<(u+75)), w, 'uni', 0), num2cell(Tmdir.SACCADE), Tmdir.(sprintf('spiketimes_%d',prb)), 'uni', 0);
            sac_fr = cell2mat(vertcat(sac_fr{:}));

            spont_fr = cellfun(@(u,w) cellfun(@(q) sum(q>=(u(1)-200) & q<(u(1)-50)), w, 'uni', 0), Tmdir.TARG_ON, Tmdir.(sprintf('spiketimes_%d',prb)), 'uni', 0);
            spont_fr = cell2mat(vertcat(spont_fr{:}));

            dp_vis = (mean(vis_fr,1) - mean(spont_fr,1)) ./ (sqrt(0.5*(var(vis_fr,1) + var(spont_fr,1))));
            dp_sac = (mean(sac_fr,1) - mean(spont_fr,1)) ./ (sqrt(0.5*(var(sac_fr,1) + var(spont_fr,1))));

            S.kilosort(prb).clusters.dp_vis = dp_vis';
            S.kilosort(prb).clusters.dp_sac = dp_sac';

            %
            dirs = sort(unique(Tmdir.angle));
            fr_perDir = zeros(numel(dirs),size(delay_fr,2));
            for d = 1:numel(dirs)
                fr_perDir(d,:) = mean(delay_fr(Tmdir.angle==dirs(d),:),1);
            end

            S.kilosort(prb).clusters.mdir_delayFR_perDir = num2cell(fr_perDir');
            [ii,mm] = max(fr_perDir);
            S.kilosort(prb).clusters.mdir_delayFR_peakDirFR = ii';
            S.kilosort(prb).clusters.mdir_delayFR_peakDir = dirs(mm);

            % VISUAL
            [vis_sel_dir, vis_pref_dir, ~, ~, frs_perAng_vis] = calculate_direction_tuning_from_tbl(Tmdir,'FR_WIN',[50,150],'ALIGN_TO','stim','PROBE',prb);
            S.kilosort(prb).clusters.vis_sel_dir = vis_sel_dir;
            S.kilosort(prb).clusters.vis_pref_dir = vis_pref_dir;

            % MOTOR
            [sac_sel_dir, sac_pref_dir, ~, ~, frs_perAng_sac] = calculate_direction_tuning_from_tbl(Tmdir,'FR_WIN',[-50,50],'ALIGN_TO','sacc','PROBE',prb);
            S.kilosort(prb).clusters.sac_sel_dir = sac_sel_dir;
            S.kilosort(prb).clusters.sac_pref_dir = sac_pref_dir;

            % VMI
            [VMI_per_unit,VMIdp_per_unit, visFR_per_unit, sacFR_per_unit] = deal(zeros(size(frs_perAng_sac,2),1));
            for unit = 1:size(frs_perAng_sac,2)
                visFR = frs_perAng_vis(:,unit);
                sacFR = frs_perAng_sac(:,unit);

                visFR_per_unit(unit) = mean(vertcat(visFR{:}));
                sacFR_per_unit(unit) = mean(vertcat(sacFR{:}));

                VMI_per_unit(unit) = (mean(vertcat(visFR{:})) - mean(vertcat(sacFR{:})))/(mean(vertcat(visFR{:})) + mean(vertcat(sacFR{:})));
                VMIdp_per_unit(unit) = (mean(vertcat(visFR{:}))-mean(vertcat(sacFR{:})))/sqrt(var(vertcat(visFR{:})) * var(vertcat(sacFR{:})));
            end
            S.kilosort(prb).clusters.vis_meanFR = visFR_per_unit;
            S.kilosort(prb).clusters.sac_meanFR = sacFR_per_unit;
            S.kilosort(prb).clusters.VMI = VMI_per_unit;
            S.kilosort(prb).clusters.VMIdp = VMIdp_per_unit;
        end
    end

    % PURSUIT
    matchingFields2 = fields(contains(fields, {'pursuit', 'purs'}, 'IgnoreCase', true));

    if ~isempty(matchingFields2)
        fprintf('\n~~CALCULATING PURS METRICS~~\n');
        Tpurs = []; 
        for mm = 1:numel(matchingFields2)
            Tpurs = [Tpurs; S.(matchingFields2{mm}).tbl];
        end

        for prb = 1:size(S.kilosort,1)
            % Percent trials fired on
            vals = cellfun(@(u) cellfun(@(q) numel(q), u, 'uni', 1)>0, Tpurs.(sprintf('spiketimes_%d',prb)), 'uni', 0);
            S.kilosort(prb).clusters.ratio_purs_trials = sum(vertcat(vals{:}),1)'./(height(Tpurs));

            % delay period FR
            targ_fr = cellfun(@(u,v,w) cellfun(@(q) (sum(q>=u & q<v)./(v-u(1)))*1000, w, 'uni', 0), num2cell(Tpurs.PURSUIT_TARG_ON), num2cell(Tpurs.PURSUIT_TARG_OFF), Tpurs.(sprintf('spiketimes_%d',prb)), 'uni', 0);
            targ_fr = cell2mat(vertcat(targ_fr{:}));

            % dprime between evoked and spotantaneous activity
            pur_fr = cellfun(@(u,w) cellfun(@(q) sum(q>=(u+50) & q<(u+250)), w, 'uni', 0), num2cell(Tpurs.PURSUIT_TARG_ON), Tpurs.(sprintf('spiketimes_%d',prb)), 'uni', 0);
            pur_fr = cell2mat(vertcat(pur_fr{:}));

            spont_fr = cellfun(@(u,w) cellfun(@(q) sum(q>=(u-200) & q<(u-50)), w, 'uni', 0), num2cell(Tpurs.PURSUIT_TARG_ON), Tpurs.(sprintf('spiketimes_%d',prb)), 'uni', 0);
            spont_fr = cell2mat(vertcat(spont_fr{:}));

            dp_pur = (mean(pur_fr,1) - mean(spont_fr,1)) ./ (sqrt(0.5*(var(pur_fr,1) + var(spont_fr,1))));

            S.kilosort(prb).clusters.dp_pur = dp_pur';

            dirs = sort(unique(Tpurs.angle));
            fr_perDir = zeros(numel(dirs),size(targ_fr,2));
            for d = 1:numel(dirs)
                fr_perDir(d,:) = mean(targ_fr(Tpurs.angle==dirs(d),:),1);
            end

            S.kilosort(prb).clusters.purs_targFR_perDir = num2cell(fr_perDir');
            [ii,mm] = max(fr_perDir);
            S.kilosort(prb).clusters.purs_targFR_peakDirFR = ii';
            S.kilosort(prb).clusters.purs_targFR_peakDir = dirs(mm);

            % MOTOR (PURSUIT)
            [pur_sel_dir, pur_pref_dir, ~, ~, ~] = calculate_direction_tuning_from_tbl(Tpurs,'FR_WIN',[-50,50],'ALIGN_TO','purs','PROBE',prb);
            S.kilosort(prb).clusters.pur_sel_dir = pur_sel_dir;
            S.kilosort(prb).clusters.pur_pref_dir = pur_pref_dir;
        end
    end

    % SPI, if dirmem and pursuit were ran
    if ~isempty(matchingFields1) & ~isempty(matchingFields2)
        for prb = 1:size(S.kilosort,1)
            % MOTOR (SACCADE)
            [~, ~, ~, ~, frs_perAng_sac] = calculate_direction_tuning_from_tbl(Tmdir,'FR_WIN',[-50,50],'ALIGN_TO','sacc','PROBE',prb);

            % MOTOR (PURSUIT)
            [~, ~, ~, ~, frs_perAng_pur] = calculate_direction_tuning_from_tbl(Tpurs,'FR_WIN',[-50,50],'ALIGN_TO','purs','PROBE',prb);

            % VMI
            [SPI_per_unit,SPIdp_per_unit,purFR_per_unit] = deal(zeros(size(frs_perAng_sac,2),1));
            for unit = 1:size(frs_perAng_sac,2)
                sacFR = frs_perAng_sac(:,unit);
                purFR = frs_perAng_pur(:,unit);

                purFR_per_unit(unit) = mean(vertcat(purFR{:}));
            
                SPI_per_unit(unit) = (mean(vertcat(sacFR{:})) - mean(vertcat(purFR{:})))/(mean(vertcat(sacFR{:})) + mean(vertcat(purFR{:})));
                SPIdp_per_unit(unit) = (mean(vertcat(sacFR{:}))-mean(vertcat(purFR{:})))/sqrt(var(vertcat(sacFR{:})) * var(vertcat(purFR{:})));
            end
            S.kilosort(prb).clusters.pur_meanFR = purFR_per_unit;
            S.kilosort(prb).clusters.SPI = SPI_per_unit;
            S.kilosort(prb).clusters.SPIdp = SPIdp_per_unit;
        end
    end

    if ~isempty(SAVE_NAME)
        save(table_path, 'S', '-v7.3');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end