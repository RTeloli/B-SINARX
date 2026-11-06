function filtered_signal =  butterWorthFilterHP(signal, dt, fc)
 
fs = 1/dt;
[b,a] = butter(9,fc/(fs/2), "low");
filtered_signal = filtfilt(b, a, signal );


end