
function [x_resampled, FsNew] = upsampler2(signal, dt, upsampling_factor)

[N, n_signal]= size(signal) ;
Fs = 1/dt;

[p, q] = rat(upsampling_factor); % Rational approximation of resampling ratio

% Resample the signal
x_resampled = resample(signal, p, q);
FsNew = upsampling_factor*Fs;
 end

