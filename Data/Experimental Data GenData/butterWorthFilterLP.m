function filtered_signal =  butterWorthFilterLP(signal, dt, fc)
 
fs = 1/dt;
[b,a] = butter(2,fc/(fs/2), "high");
filtered_signal = filtfilt(b, a, signal );


end