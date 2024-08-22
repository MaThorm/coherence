% Script that takes the xx_trials script and exchanges them for artificial
% data

function [new_trials] = do_create_data(trials,freqs,amp_range,noise_f)
% amp range(1) = lower amp bound
% amp range(2) = upper amp bound
% noise f = noise amplifcation factor
new_trials = trials;
Fs = 1000;


for i_t = 1:length(trials)
    for i_s = 1:length(trials(i_t).trial)
        for i_c = 1:size(trials(i_t).trial{i_s},1)
            cur_length = length(trials(i_t).trial{i_s}(i_c,:));
            signal = zeros(1,cur_length);
            t = 1:cur_length;
            for i_f = 1:length(freqs)
                signal = signal + sin(2*pi*freqs(i_f)*t);
            end 
            amp_mod = (rand(1) * (amp_range(2)-amp_range(1)) + amp_range(1));
            signal = signal * amp_mod;
            noise = randn(1,cur_length) * noise_f;
            signal = signal + noise;
            new_trials(i_t).trial{i_s}(i_c,:) = signal;
        end 
    end 
end 
end 
