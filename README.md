# Sample_Entropy

"SampEn compares segments of the time series to a template of length m + 1. If the first m timepoints match the template (within a tolerance factor r) the segment is listed as an “m match.” If all m + 1 timepoints match the template within the tolerance then the segment is also listed as an “m + 1 match.” The template-matching process is repeated so that each segment is considered a template once, and is also assessed for matching the other segments many times. The proportion of m + 1 matches to m matches is considered a measure of complexity (i.e., if a high proportion of the length m matches are also length m + 1 matches, then the time series is predictable and has low complexity). SampEn is the negative log of this proportion. As described by Lake et al. (6), higher m values and lower r values tend to reduce both the number of length m matches (Cm) and the number of length m + 1 matches (Cm1)." <br/>

Roediger, D. J., Butts, J., Falke, C., Fiecas, M. B., Klimes-Dougan, B., Mueller, B. A., & Cullen, K. R. (2024). Optimizing the measurement of sample entropy in resting-state fMRI data. Frontiers in Neurology, 15, 1331365.

![image](https://github.com/user-attachments/assets/2e46c21a-3ea8-4b2a-961d-dc56fe156e3c)

## Code (ChatGPT)

The implementation closely follows the method described by Richman and Moorman in their original paper on Sample Entropy:
Richman, J.S., & Moorman, J.R. (2000). "Physiological time-series analysis using approximate entropy and sample entropy." American Journal of Physiology-Heart and Circulatory Physiology, 278(6), H2039-H2049.

### (MATLAB)

    function sampen_value = SampEn(data, m, r)
    
        % data: time series data (1D array)
        % m: embedding dimension
        % r: tolerance (often set as r = 0.2 * std(data))
    
        N = length(data);  % Length of the time series
        % Step 1: Create vectors of length m
        X_m = zeros(N - m + 1, m);
        for i = 1:(N - m + 1)
            X_m(i, :) = data(i:i + m - 1);
        end
        
        % Step 2: Count the number of matches within the tolerance r
        B = 0;  % Similarity count for length m
        A = 0;  % Similarity count for length m+1
    
        for i = 1:(N - m + 1)
            % Distance calculation for vectors of length m
            distance = max(abs(X_m(i+1:end, :) - X_m(i, :)), [], 2);
            B = B + sum(distance < r);
        end
    
        % Create vectors of length m+1
        X_m1 = zeros(N - m, m + 1);
        for i = 1:(N - m)
            X_m1(i, :) = data(i:i + m);
        end
        
        for i = 1:(N - m)
            % Distance calculation for vectors of length m+1
            distance = max(abs(X_m1(i+1:end, :) - X_m1(i, :)), [], 2);
            A = A + sum(distance < r);
        end
    
        % Step 3: Calculate SampEn
        B = B / (N - m + 1);
        A = A / (N - m);
    
        % SampEn is the negative natural logarithm of the ratio of A to B
        sampen_value = -log(A / B);
        
    end

#### Example of use:

    % Example time series data
    data = randn(1, 1000); % Random data with 1000 points
    
    % Parameters
    m = 2;
    r = 0.2 * std(data);
    
    % Calculate SampEn
    sampen_value = SampEn(data, m, r);
    disp(['Sample Entropy: ', num2str(sampen_value)]);

#### Explanation:
Embedding Dimension (m): This is the length of the sequences to be compared. Typically, m is set to 2. </br>

Tolerance (r): This is a threshold distance for considering two sequences as similar. A common choice is r = 0.2 * std(data). </br>

Step 1: Forming Vectors: </br>
The data is divided into overlapping vectors of length m.</br>

Step 2: Counting Matches: </br>
For each vector, count the number of other vectors within the distance r.</br>
This is done first for vectors of length m and then for vectors of length m + 1.</br>

Step 3: Calculating Sample Entropy: </br>
The ratio of counts for length m + 1 and m is calculated.</br>
SampEn is the negative logarithm of this ratio.</br>

### Python

    import numpy as np
    
    def sample_entropy(data, m, r):
        """
        Calculate the Sample Entropy (SampEn) of a time series.
    
        Parameters:
        data : list or numpy array
            The time series data (1D array).
        m : int
            The embedding dimension (length of sequences to be compared).
        r : float
            The tolerance (usually set as r = 0.2 * std(data)).
    
        Returns:
        sampen_value : float
            The calculated Sample Entropy value.
        """
        N = len(data)  # Length of the time series
    
        # Step 1: Create vectors of length m
        X_m = np.array([data[i:i+m] for i in range(N - m + 1)])
    
        # Step 2: Count the number of matches within the tolerance r
        B = 0  # Similarity count for length m
        A = 0  # Similarity count for length m+1
    
        for i in range(N - m + 1):
            # Distance calculation for vectors of length m
            distance = np.max(np.abs(X_m[i+1:] - X_m[i]), axis=1)
            B += np.sum(distance < r)
    
        # Create vectors of length m+1
        X_m1 = np.array([data[i:i+m+1] for i in range(N - m)])
    
        for i in range(N - m):
            # Distance calculation for vectors of length m+1
            distance = np.max(np.abs(X_m1[i+1:] - X_m1[i]), axis=1)
            A += np.sum(distance < r)
    
        # Step 3: Calculate SampEn
        B = B / (N - m + 1)
        A = A / (N - m)
    
        # SampEn is the negative natural logarithm of the ratio of A to B
        sampen_value = -np.log(A / B)
    
        return sampen_value
    
    # Example of Use
    data = np.random.randn(1000)  # Random data with 1000 points
    
    # Parameters
    m = 2
    r = 0.2 * np.std(data)
    
    # Calculate SampEn
    sampen_value = sample_entropy(data, m, r)
    print(f"Sample Entropy: {sampen_value:.4f}")


## Toolboxes

# Open questions:
1.- If r is calculated on each segment separately, then the increase of non-EEG noise (for example) could alter the measurement of SampEn. Should we calulate r on the entire experiment instead of on each trial? I would do it separately for each electrode.




