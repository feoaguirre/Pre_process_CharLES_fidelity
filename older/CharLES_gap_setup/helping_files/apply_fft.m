% Given time and signal vectors of equal even length, this function applies fft and returns the
% peak amplitudes and their frequencies (regular frequencies in 1/[time], not angular frequency). 
% threshhold_value allows you to get rid of peaks in the spectrum created by noise.

function fft_results = apply_fft(time, signal, threshhold_value)
    % Parameters
    L = length(signal); % Number of samples
    Fs = L/(time(end) - time(1)); % Sampling frequency 1/[time]

    % Getting FFT results
    pure_fft = fft(signal);
    normalized_fft = abs(pure_fft/L);
    single_sided_fft = normalized_fft(1:L/2+1);
    single_sided_fft(2:end-1) = 2*single_sided_fft(2:end-1);
    single_sided_freqs = Fs/L*(0:(L/2));
    
    % Getting rid of noise
    single_sided_fft(single_sided_fft < threshhold_value) = 0;
    
    % Finding peaks
    peak_indices = find(islocalmax(single_sided_fft));
    if single_sided_fft(1) > 0
        peak_indices = [1; peak_indices];
    end
    fft_results.frequencies = single_sided_freqs; % Units of 1/[time]
    fft_results.amplitudes = single_sided_fft;
    fft_results.peak_frequencies = single_sided_freqs(peak_indices);
    fft_results.peak_amplitudes = single_sided_fft(peak_indices);
end
