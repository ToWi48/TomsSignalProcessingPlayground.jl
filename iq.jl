using FFTW
using PlotlyJS

Fs = 64
T = 1/Fs
N = 512
t = (0:N-1)*T

f = -1
phi = deg2rad(0)
f0 = 2
phi0 = deg2rad(0)

s = exp.(im*(2*pi*f*t.+phi))
s_spec = fftshift(1 / length(s) * fft(s))

display(plot(
    [
        scatter3d(x=t, y=real(s), z = imag(s), name="Signal (Imag)", mode="lines"), 
    ], 
    Layout(
        title="Zeitsignal",
        scene = attr(
            xaxis = attr(title = "Time [s]"),
            yaxis = attr(title = "Realpart", range=[-1, 1]),
            zaxis = attr(title = "Imaginarypart", range=[-1, 1])
        )
    )
))

display(plot(
    [
        scatter(x=Fs/N*(-N/2:N/2-1), y=abs.(s_spec), name="Signal"), 
    ], 
    Layout(
        title="Zeitsignal",
        xaxis_title="Frequency (Hz)",
        yaxis_title="Amplitude",
        yaxis=attr(range=[0, 1])
    )
))

# Modulator
si = cos.(2*pi*f0*t.+phi0) .* real(s)
sq = -sin.(2*pi*f0*t.+phi0) .* imag(s)

sm = si + sq
sm_spec = fftshift(1 / length(sm) * fft(sm))

display(plot(
    [
        scatter3d(x=t, y=real(sm), z = imag(sm), name="Signal (Imag)", mode="lines"), 
    ], 
    Layout(
        title="Modulatorausgangssignal",
        scene = attr(
            xaxis = attr(title = "Frequency (Hz)"),
            yaxis = attr(title = "Realpart", range=[-1, 1]),
            zaxis = attr(title = "Imaginarypart", range=[-1, 1])
        )
    )
))

display(plot(
    [
        scatter(x=Fs/N*(-N/2:N/2-1), y=abs.(sm_spec), name="Signal"), 
    ], 
    Layout(
        title="Modulatorausgangssignal Spec",
        xaxis_title="Frequency (Hz)",
        yaxis_title="Amplitude",
        yaxis=attr(range=[0, 1])
    )
))

# Demodulator
sdmi = cos.(2*pi*f0*t) .* sm
sdmq = -sin.(2*pi*f0*t) .* sm

sdmik = sdmi .- 0.5*cos.(2*pi*(2*f0+f)*t.+phi0.+phi)
sdmqk = sdmq .+ 0.5*sin.(2*pi*(2*f0+f)*t.+phi0.+phi)

sdm = sdmik + im * sdmqk
sdm_spec = fftshift(1 / length(sdm) * fft(sdm))

display(plot(
    [
        scatter3d(x=x=Fs/N*(-N/2:N/2-1), y=1/N*real(fftshift(fft(sdmi))), z = zeros(N), name="Inphase", mode="lines"), 
        scatter3d(x=x=Fs/N*(-N/2:N/2-1), y=zeros(N), z=1/N*imag(fftshift(fft(sdmq))), name="Quadratur", mode="lines"), 
    ], 
    Layout(
        title="Demodulatorausgangssignal",
        scene = attr(
            xaxis = attr(title = "Frequency (Hz)"),
            yaxis = attr(title = "Realpart", range=[-1, 1]),
            zaxis = attr(title = "Imaginarypart", range=[-1, 1])
        )
    )
))

display(plot(
    [
        scatter(x=Fs/N*(-N/2:N/2-1), y=abs.(sdm_spec), name="Signal"), 
    ], 
    Layout(
        title="Demodulatorausgangssignal Spec",
        xaxis_title="Frequency (Hz)",
        yaxis_title="Amplitude",
        yaxis=attr(range=[0, 1])
    )
))