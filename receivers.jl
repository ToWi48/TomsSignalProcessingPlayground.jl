using FFTW
using PlotlyJS

function Lowpass(s, fg, Fs)
    S = fft(s)                          # transform signal from time to frequency domain
    S[Int(end/2):end] .= 0              # remove imaginary part
    start_bin = Int(length(s)*fg/Fs)    # calc bin index depending on cutoff frequency
    S[start_bin:end] .= 0               # filter by cutoff index
    return 2*real(ifft(S))              # transform signal back from frequency to time domain
end

function Highpass(s, fg, Fs)
    S = fft(s)                          # transform signal from time to frequency domain
    S[Int(end/2):end] .= 0              # remove imaginary part
    stop_bin = Int(length(s)*fg/Fs)     # calc bin index depending on cutoff frequency
    S[1:stop_bin] .= 0                  # filter by cutoff index
    return 2*real(ifft(S))              # transform signal back from frequency to time domain
end

function Hilbert(s)                     # make real signal complex
    S = 2*fft(s)                        # transform real signal to frequency domain
    S[Int(end/2):end] .= 0              # remove negative frequencies
    return ifft(S)                      # transform back into time domain
end

# Setup
Fs = 64
T = 1/Fs
N = 512
t = (0:N-1)*T

p_lo = 0
f_lo = 10
s_lo = cos.(2*pi*f_lo*t.+deg2rad(p_lo))
s_lo_spec = fftshift(1 / length(s_lo) * fft(s_lo))

f_if = 3

p_rf = 0
f_rf = f_lo+f_if
s_rf = cos.(2*pi*f_rf*t.+deg2rad(p_rf))
s_rf_spec = fftshift(1 / length(s_rf) * fft(s_rf))

p_rf_img = 0
f_rf_img = f_lo-f_if+1
s_rf_img = cos.(2*pi*f_rf_img*t.+deg2rad(p_rf_img))
s_rf_img_spec = fftshift(1 / length(s_rf_img) * fft(s_rf_img))

s_in = s_rf #+ s_rf_img

# plot results
plot_time = plot(
    [
        scatter(x=Fs/N*(-N/2:N/2-1), y=abs.(s_lo_spec), name="Local Oscillator"), 
        scatter(x=Fs/N*(-N/2:N/2-1), y=abs.(s_rf_spec), name="RF Signal"), 
        scatter(x=Fs/N*(-N/2:N/2-1), y=abs.(s_rf_img_spec), name="RF* Signal")
    ], 
    Layout(
        title="Source Signals",
        xaxis_title="Frequency (Hz)",
        yaxis_title="Amplitude",
        yaxis=attr(range=[0, 1])
    )
)

########################################################
# Heterodyne Receiver
########################################################

s_mixed = s_in .* s_lo
s_mixed_spec = fftshift(1 / length(s_mixed) * fft(s_mixed))

# Remix signal, but with image rejection filter
s_mixed_image_reject = Highpass(s_in, f_lo, Fs) .* s_lo

# Lowpass filter to filter out the addition term
s_mixed_lowpass = Lowpass(s_mixed, f_lo+f_if, Fs)
s_mixed_image_reject_lowpass = Lowpass(s_mixed_image_reject, f_lo+f_if, Fs)

s_mixed_lowpass_spec = fftshift(1 / length(s_mixed_lowpass) * fft(s_mixed_lowpass))
s_mixed_image_reject_lowpass_spec = fftshift(1 / length(s_mixed_image_reject_lowpass) * fft(s_mixed_image_reject_lowpass))

# channel selection skiped -> same procedure as every year, james!

# plot results
plot_freq = plot(
    [
        scatter(x=Fs/N*(-N/2:N/2-1), y=abs.(s_mixed_spec), name="Mixer output"), 
        scatter(x=Fs/N*(-N/2:N/2-1), y=abs.(s_mixed_lowpass_spec), name="Mixer output Lowpass"),
        scatter(x=Fs/N*(-N/2:N/2-1), y=abs.(s_mixed_image_reject_lowpass_spec), name="Mixer output Lowpass (Image Reject)"),
    ],
    Layout(
        title="Heterodyne Receiver",
        xaxis_title="Frequency (Hz)",
        yaxis_title="Amplitude",
        yaxis=attr(range=[0, 1])
    )
)

########################################################
# Hartley Receiver
########################################################

### upper path
s_lo_sin = real(Hilbert(s_lo) .* exp(im * deg2rad(-90)))  # somehow this receiver uses a sin instead of a cosin oscillator, so simply return the cosine LO back by 90 degrees
s_hartley_upper_mixed = s_in .* s_lo_sin  
s_hartley_upper_mixed_filtered = Lowpass(s_hartley_upper_mixed, f_lo+f_if, Fs)
s_hartley_upper = real(Hilbert(s_hartley_upper_mixed_filtered) .* exp(im * deg2rad(-90)))
s_hartley_upper_spec = fftshift(1 / length(s_hartley_upper) * fft(s_hartley_upper))

### lower path
s_hartley_lower_mixed = s_in .* s_lo    # because s_lo is real cosine
s_hartley_lower = Lowpass(s_hartley_lower_mixed, f_lo+f_if, Fs)
s_hartley_lower_spec = fftshift(1 / length(s_hartley_lower) * fft(s_hartley_lower))

# merge both paths together
s_hartley_out = s_hartley_upper .+ s_hartley_lower
s_hartley_out_spec = fftshift(1 / length(s_hartley_out) * fft(s_hartley_out))

# plot results
plot_hartley_freq = plot(
    [
        scatter(x=Fs/N*(-N/2:N/2-1), y=abs.(s_hartley_upper_spec), name="Hartley Upper Path"), 
        scatter(x=Fs/N*(-N/2:N/2-1), y=abs.(s_hartley_lower_spec), name="Hartley Lower Path"),
        scatter(x=Fs/N*(-N/2:N/2-1), y=abs.(s_hartley_out_spec), name="Hartley Out"),
    ],
    Layout(
        title="Hartley Receiver",
        xaxis_title="Frequency (Hz)",
        yaxis_title="Amplitude",
        yaxis=attr(range=[0, 1])
    )
)

########################################################
# Direct Conversion
########################################################

s_direct_out = s_in .* s_rf   # s_lo = s_rf
s_direct_out_spec = fftshift(1 / length(s_direct_out) * fft(s_direct_out))

# plot results
plot_direct_freq = plot(
    [
        scatter(x=Fs/N*(-N/2:N/2-1), y=abs.(s_direct_out_spec), name="Direct Conversion out"), 
    ],
    Layout(
        title="Direct Conversion",
        xaxis_title="Frequency (Hz)",
        yaxis_title="Amplitude",
        yaxis=attr(range=[0, 1])
    )
)

########################################################
# IQ (analog, digital (QAM) is boring...)
########################################################

# Inphase path
s_i_mixed = 2 * s_in .* s_lo    # s_lo = cosine
s_i_mixed = Lowpass(s_i_mixed, f_rf+1, Fs)
s_i_mixed_spec = fftshift(1 / length(s_i_mixed) * fft(s_i_mixed))

# Quadratur path
s_q_mixed = 2 * s_in .* real(Hilbert(s_lo) .* exp(im * deg2rad(90)))    # s_lo = -sin
s_q_mixed = Lowpass(s_q_mixed, f_rf+1, Fs)
s_q_mixed_spec = fftshift(1 / length(s_q_mixed) * fft(s_q_mixed))

s_iq_out = s_i_mixed + s_q_mixed * im
s_iq_out_spec = fftshift(1 / length(s_iq_out) * fft(s_iq_out))

# plot results
plot_iq_freq = plot(
    [
        scatter3d(x=x=Fs/N*(-N/2:N/2-1), y=1/N*real(fftshift(fft(s_i_mixed))), z = zeros(N), name="Inphase", mode="lines"), 
        scatter3d(x=x=Fs/N*(-N/2:N/2-1), y=zeros(N), z=1/N*imag(fftshift(fft(s_q_mixed))), name="Quadratur", mode="lines"),
    ],
    Layout(
        title="IQ Demodulator",
        xaxis=attr(title="Frequency (Hz)"),
        yaxis=attr(title="Real Part", range=[0, 1]),
        zaxis=attr(title="Imaginary Part", range=[0, 1]),
        scene=attr(
            xaxis_title="Real Part",
            yaxis_title="Imaginary Part",
            zaxis_title="Frequency (Hz)"
        )
    )
)

plot_iq_spec = plot(
    [
        scatter(x=Fs/N*(-N/2:N/2-1), y=abs.(s_iq_out_spec), name="Signal"), 
    ], 
    Layout(
        title="Demodulatorausgangssignal Spec",
        xaxis_title="Frequency (Hz)",
        yaxis_title="Amplitude",
        yaxis=attr(range=[0, 1])
    )
)

########################################################
# Plotting
########################################################

fig = [
    plot_time plot_freq plot_hartley_freq plot_direct_freq;
    plot_iq_freq plot_iq_spec
]
relayout!(fig)

display(fig)