using PlotlyJS
using FFTW

# Zeitsignal
Fs = 64
T = 1/Fs
N = 512
t = (0:N-1)*T

p_lo = 0
f_lo = 1
s_lo = exp.(im*2*pi*f_lo*t.+deg2rad(p_lo))
s_lo_spec = fftshift(1 / length(s_lo) * fft(s_lo))

plot_time = plot(
    [
        scatter(x=t, y=real(s_lo), name="Signal (Real)"), 
        scatter(x=t, y=imag(s_lo), name="Signal (Imag)"), 
    ], 
    Layout(
        title="Zeitsignal",
        xaxis_title="Frequency (Hz)",
        yaxis_title="Amplitude",
        yaxis=attr(range=[-1, 1])
    )
)

plot_freq = plot(
    [
        scatter(x=Fs/N*(-N/2:N/2-1), y=abs.(s_lo_spec), name="Signal"), 
    ], 
    Layout(
        title="Zeitsignal",
        xaxis_title="Frequency (Hz)",
        yaxis_title="Amplitude",
        yaxis=attr(range=[0, 1])
    )
)

fig = [
    plot_time plot_freq;
]
relayout!(fig)

display(fig)