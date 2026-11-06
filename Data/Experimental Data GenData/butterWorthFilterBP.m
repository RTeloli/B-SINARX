function filtered_signal =  butterWorthFilterBP(signal, dt, fc)
wn = [15 fc];
fs = 1/dt;
[b,a] = butter(9,wn/(fs/2), "bandpass");
filtered_signal = filtfilt(b, a, signal );


end