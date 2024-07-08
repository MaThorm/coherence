function [new_trials,inc,test] = do_SSD(trials,fs,th,lower,upper)
%DO_SSD Function that takes a trialselected fieldtrip file and, using that,
%switches out each trial with the biggest component of the SSD function.
% trials: fieldtrip result of do_trialselection
% Fs: sampling frequency
% th: threshold which defines the variance of the residual (default
% 0.01),
% Upper and lower frequency boundary
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Gives back
% new_trials: the new trial block 
% inc: whether the trial, according to the set criteria will be included or
% not
new_trials = trials;
parfor i_sess = 1:length(trials)
    for i_t = 1:length(trials(i_sess).trial)
        for i_c = 1:size(trials(i_sess).trial{1, i_t},1)
            cur_trial = trials(i_sess).trial{1, i_t}(i_c,:);
            cur_SSD = SSD(cur_trial,fs,th,10); %does the SSD algorithm on a single trial MAX 10 components right now
            test(i_sess).SSD{i_t,i_c} = cur_SSD;
%             [p,f] = pspectrum(cur_SSD',fs); % calculates the powerspectrum of each SSD component
            
            [p, ntaper, f] = ft_specest_mtmfft(cur_SSD,0.001:0.001:length(cur_trial)/1000,'taper','hanning');
            % P seems to be the c x f in the case where there is only one
            % component, adding a condition which changes that 
            if size(cur_SSD,1) > 1            
                p = abs(squeeze(p))';
            else 
                p = abs(squeeze(p));
            end 
            f = f';
            test(i_sess).p{i_t,i_c} = p;
            test(i_sess).f{i_t,i_c} = f';
            l_b = min(find(f >= lower)); % Finds the frequency value of lower bound
            u_b = max(find(f <= upper)); % Same for upper
            full = trapz(p);
            test(i_sess).full{i_t,i_c} = full;
            bounded = trapz(p(l_b:u_b,:));     % Calculates area under curve for all powerspectra within upper and lower bound
            test(i_sess).bounded{i_t,i_c} = bounded;
            q = bounded./full;
            test(i_sess).q{i_t,i_c} = q;
            %[a,m_ex] = max(q)
            [a,m_ex] = max(bounded);    
            new_trials(i_sess).trial{1, i_t}(i_c,:) = cur_SSD(m_ex,:);
            [a,peak] = max(p(:,m_ex));
            peak_freq = f(peak);
            % Testing three things: whether proportion is bigger than 0.7
            % And whether peak is within frequency bounds 
            if q(m_ex) < 0.7 || peak_freq < lower || peak_freq > upper 
                inc(i_sess).inc{i_t,i_c} = false;
            else 
                inc(i_sess).inc{i_t,i_c} = true;
            end
            fprintf(sprintf("Session %d, trial %d, component %d",i_sess,i_t,i_c))
        end 
    end 
end  
end
