# using DSP
using Plots

let 

    function f(x, y)
        x[1] = 42    # mutates x
        y += 7    # new binding for y, no mutation
        return y
    end

    a = [1,2,3]
    b = 2

    f(a,b)

    b

    # # Generate a noisy signal
    # sample_rate = 1000.0  # Sample rate in Hz
    # t = 0:1/sample_rate:1  # Time vector
    # x = sin.(2π * 100 * t) + 0.5 * randn(length(t))  # Noisy signal (sinusoid + noise)

    # # Plot the original signal
    # pl = plot(t, x, label="Noisy Signal")
    # # plot!(t, filtered_signal, label="Filtered Signal")
    # xlabel!("Time (s)")
    # ylabel!("Amplitude")
    # # title!("Low-Pass Filtered Signal")
    
    # # filter
    # f_cutoff = 0.2
    # responsetype = Lowpass(10; fs=sample_rate)
    # designmethod = Butterworth(4)
    # xf = filt(digitalfilter(responsetype, designmethod), x)  

    # # Plot the filtered signal
    # plot!(t, xf, label="Filtered Signal")
    # title!("Low-Pass Filtered Signal")
    # display(pl)

    
end
