using  Plots
using Statistics

let 

    # Example data (replace this with your actual sample)
    data = [1.2, 2.5, 3.1, 4.0, 5.2]

    # Calculate the mean
    mean_value = mean(data)

    # Calculate the variance with known mean
    sample_variance_with_mean = var(data)

    # Print the result
    println("Sample Variance (given mean): ", sample_variance_with_mean)

end