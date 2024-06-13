function [new_inst,cut,testing] = do_snippeting(wavelet_SSD,inst,yrange,toi,sd_mult,outlier_mult)
%DO_SNIPPETING Based on SSD wavelet power checks if power is high enough
%for hilber transform. Gives back inst. freq. with reduced power removed
%as well as the logical array the removal is based on. 
%~~~~~~~~~~~~~~~Input~~~~~~~~~~~~~~~~~%
% wavelet_SSD = Wavelet of the SSD component with .in, .out and .V4 as
% fields
% inst = inst. freq of the SSD component with .in .out and .V4 as fields
% yrange = Gamma range used to define which parts of the wavelet further
% toi = time range
% analysis is based on 
% sd_mult = multiplier of the SD of all wavelet power values, higher sd_mult
% means less parts are removed
% outlier_mult = multiplier of the IQR of the wavelet power values, 

new_inst = inst;
for i_s = 1:length(inst.in)
    wave_in = wavelet_SSD.in(i_s);
    wave_out = wavelet_SSD.out(i_s);
    wave_V4 = wavelet_SSD.V4(i_s);
    inst_in = inst.in(i_s);
    inst_out = inst.out(i_s);
    inst_V4 = inst.V4(i_s);
    
    clear pow_in pow_out pow_V4
    % Gettingn only the relevant 
    for i_t = 1:size(wave_in.powspctrm,1)
        pow_in(i_t,:,:) = squeeze(wave_in.powspctrm(i_t,:,yrange,toi(1):toi(2)));
%         pow_in(i_t,:,:) = squeeze(wave_in.powspctrm(i_t,:,:,toi(1):toi(2)));
    end 
    
    for i_t = 1:size(wave_out.powspctrm,1)
        pow_out(i_t,:,:) = squeeze(wave_out.powspctrm(i_t,:,yrange,toi(1):toi(2)));
%         pow_out(i_t,:,:) = squeeze(wave_out.powspctrm(i_t,:,:,toi(1):toi(2)));
    end 
    
    for i_t = 1:size(wave_V4.powspctrm,1)
        pow_V4(i_t,:,:) = squeeze(wave_V4.powspctrm(i_t,:,yrange,toi(1):toi(2)));
%         pow_V4(i_t,:,:) = squeeze(wave_V4.powspctrm(i_t,:,:,toi(1):toi(2)));
    end 
    
    % Reshaping the gamma wavelet power spectra to 1D
    pow_1d_in = reshape(pow_in,1,[]); 
    pow_1d_in = pow_1d_in(~isnan(pow_1d_in));
      
    pow_1d_out = reshape(pow_in,1,[]); 
    pow_1d_out = pow_1d_out(~isnan(pow_1d_out));
    
    pow_1d_comb = [pow_1d_in pow_1d_out];
    testing.V1(i_s).full_pow = pow_1d_comb;

    pow_1d_V4 = reshape(pow_V4,1,[]);
    pow_1d_V4 = pow_1d_V4(~isnan(pow_1d_in));
    testing.V4(i_s).full_pow = pow_1d_V4;

    % Statistics old
    % p_m = mean(pow_1d_comb);
    % p_sd = std(pow_1d_comb);
    
    % Creating a cutoff. To reduce influence of outliers to the mean and sd,
    % the outliers larger than the median + 1.5*the iqr are cut out. This is done for each
    % electrode seperately. Meaning that both V1a and V1n get pooled together.
    % For V1
    p_med = median(pow_1d_comb,'omitnan');
    pow_iqr = iqr(pow_1d_comb);
    testing.V1(i_s).out_cut = p_med + outlier_mult*pow_iqr;
    pow_1d_comb = pow_1d_comb(pow_1d_comb <= p_med + outlier_mult*pow_iqr);
    p_m = mean(pow_1d_comb,'omitnan');
    p_sd = std(pow_1d_comb,'omitnan');
    % Cutoff is then the mean - X * the SD
    cut_off_comb = p_m - sd_mult*p_sd;
    testing.V1(i_s).pow = pow_1d_comb;
    testing.V1(i_s).p_mean = p_m;
    testing.V1(i_s).p_sd = p_sd;
    testing.V1(i_s).cut_off = cut_off_comb;
    
    % For V4
    p_med_V4 = median(pow_1d_V4,'omitnan');
    pow_iqr_V4 = iqr(pow_1d_V4);
    testing.V4(i_s).out_cut = p_med_V4 + outlier_mult*pow_iqr_V4;
    pow_1d_V4 = pow_1d_V4(pow_1d_V4 <= p_med_V4 + outlier_mult*pow_iqr_V4);
    p_mV4 = mean(pow_1d_V4,'omitnan');
    p_sdV4 = std(pow_1d_V4,'omitnan');
    cut_off_V4 = p_mV4 - sd_mult*p_sdV4;
    testing.V4(i_s).pow = pow_1d_V4;
    testing.V4(i_s).p_mean = p_mV4;
    testing.V4(i_s).p_sd = p_sdV4;
    testing.V4(i_s).cut_off = cut_off_V4;
    
    % Creating logical arrays (not actually logical) for each gamma power
    % spectra. For each time point it is evaluated whether the peak power in
    % the gamma range is smaller than the cutoff. If yes, the cut array will be
    % NaN at that point, otherwise it will be 1. Also if no peak power
    % exists, the result is also nan. 
    for i_t = 1:size(pow_in,1)
        for i_time = 1:size(pow_in,3)
            if max(pow_in(i_t,:,i_time)) <= cut_off_comb
                cut_in(i_s,i_t,i_time) = nan;
            elseif isnan(pow_in(i_t,:,i_time))
                cut_in(i_s,i_t,i_time) = nan;
            else 
                cut_in(i_s,i_t,i_time) = 1;            
            end 
        end 
        cur_trial = inst.in(i_s).trial{:,i_t};
        cur_cut = squeeze(cut_in(i_s,i_t,1:length(cur_trial)));
        cur_trial(:,isnan(cur_cut)) = nan;
        new_inst.in(i_s).trial{:,i_t} = cur_trial;
    end 
    for i_t = 1:size(pow_out,1)
        for i_time = 1:size(pow_out,3)
            if max(pow_out(i_t,:,i_time)) <= cut_off_comb
                cut_out(i_s,i_t,i_time) = nan;
            elseif isnan(pow_out(i_t,:,i_time))
                cut_out(i_s,i_t,i_time) = nan;
            else 
                cut_out(i_s,i_t,i_time) = 1;            
            end 
        end 
        cur_trial = inst.out(i_s).trial{:,i_t};
        cur_cut = squeeze(cut_out(i_s,i_t,1:length(cur_trial)));
        cur_trial(:,isnan(cur_cut)) = nan;
        new_inst.out(i_s).trial{:,i_t} = cur_trial;
    end 
    
    for i_t = 1:size(pow_V4,1)
        for i_time = 1:size(pow_V4,3)
            if max(pow_V4(i_t,:,i_time)) <= cut_off_V4
                cut_V4(i_s,i_t,i_time) = nan;
            elseif isnan(pow_V4(i_t,:,i_time))
                cut_V4(i_s,i_t,i_time) = nan;
            else 
                cut_V4(i_s,i_t,i_time) = 1;            
            end 
        end 
        cur_trial = inst.V4(i_s).trial{:,i_t};
        cur_cut = squeeze(cut_V4(i_s,i_t,1:length(cur_trial)));
        cur_trial(isnan(cur_cut)) = nan;
        new_inst.V4(i_s).trial{:,i_t} = cur_trial;
    end

end 
cut.in = cut_in;
cut.out = cut_out;
cut.V4 = cut_V4;
end

