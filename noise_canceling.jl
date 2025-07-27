using FFTW, Plots, Random

function Lowpass(s, fg, Fs)
    S = fft(s)                          # transform signal from time to frequency domain
    S[Int(end/2):end] .= 0              # remove imaginary part
    start_bin = Int(length(s)*fg/Fs)    # calc bin index depending on cutoff frequency
    S[start_bin:end] .= 0               # filter by cutoff index
    return 2*real(ifft(S))              # transform signal back from frequency to time domain
end

# Parameter
fs = 10_000           # Abtastrate in Hz
T = 1.0               # Dauer in Sekunden
N = Int(fs * T)       # Anzahl der Samples
t = range(0, stop=T, length=N)


function generate_pink_noise(N::Int)
    # Weißes Rauschen erzeugen
    white_noise = randn(N)

    # Fourier-Transformation
    spectrum = fft(white_noise)

    # Frequenzskalierung: 1/f
    freqs = [max(1, i) for i in 1:N]
    scaling = 1.0 ./ sqrt.(freqs)

    # Skaliertes Spektrum
    pink_spectrum = spectrum .* scaling

    # Rücktransformation
    pink_noise = real(ifft(pink_spectrum))

    # Normalisierung auf Amplitude 1
    pink_noise ./= maximum(abs, pink_noise)

    return pink_noise
end

noise = generate_pink_noise(N)

# Chopper-Signal (Rechteckwelle)
chopper_freq = 50  # Hz
square_wave = sign.(sin.(2π * chopper_freq .* t))

# Vin
Vin = 1 * 0.1
 
# Modulation und Demodulation
modulated = Vin .* square_wave + (noise .* 0.1)
modulated_amplified = 10 .* modulated
demodulated = modulated_amplified .* square_wave

# PSD-Funktion
function plot_psd(signal, fs, title_str)
    p = plot(title=title_str, xlabel="Frequenz (Hz)", ylabel="Leistung (dB)", legend=false)
    psd = abs.(fft(signal)).^2
    freqs = (0:N-1) .* (fs / N)
    plot!(p, freqs[1:div(N,2)], 10 .* log10.(psd[1:div(N,2)] .+ 1e-12))
    return p
end

p1 = plot([modulated])
p2 = plot_psd(modulated, fs, "Nach Chopper-Modulation")
p3 = plot(Lowpass(demodulated, 3, fs))
#p2 = plot([modulated, modulated_noise])
#p1 = plot_psd(noise, fs, "Originales 1/f-Rauschen")
#p2 = plot_psd(modulated, fs, "Nach Chopper-Modulation")
#p3 = plot_psd(demodulated, fs, "Nach Demodulation")

plot(p1, p2, p3, layout=(3,1), size=(900, 800))
