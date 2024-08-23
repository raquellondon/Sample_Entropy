# Sample_Entropy

## Explanation
<details><summary>click here for more</summary><p>

"SampEn compares segments of the time series to a template of length m + 1. If the first m timepoints match the template (within a tolerance factor r) the segment is listed as an “m match.” If all m + 1 timepoints match the template within the tolerance then the segment is also listed as an “m + 1 match.” The template-matching process is repeated so that each segment is considered a template once, and is also assessed for matching the other segments many times. The proportion of m + 1 matches to m matches is considered a measure of complexity (i.e., if a high proportion of the length m matches are also length m + 1 matches, then the time series is predictable and has low complexity). SampEn is the negative log of this proportion. As described by Lake et al. (6), higher m values and lower r values tend to reduce both the number of length m matches (Cm) and the number of length m + 1 matches (Cm1)." <br/>

Roediger, D. J., Butts, J., Falke, C., Fiecas, M. B., Klimes-Dougan, B., Mueller, B. A., & Cullen, K. R. (2024). Optimizing the measurement of sample entropy in resting-state fMRI data. Frontiers in Neurology, 15, 1331365.

![image](https://github.com/user-attachments/assets/2e46c21a-3ea8-4b2a-961d-dc56fe156e3c)

</details>

## Code
<details><summary>click here for more</summary><p>

### SampEn functions
<details><summary>click here for more</summary><p>

The implementation closely follows the method described by Richman and Moorman in their original paper on Sample Entropy:
Richman, J.S., & Moorman, J.R. (2000). "Physiological time-series analysis using approximate entropy and sample entropy." American Journal of Physiology-Heart and Circulatory Physiology, 278(6), H2039-H2049.

### MATLAB - This code has been tested

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

### Python - This code has not been tested

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
</details>

### Analysis code (MATLAB)
<details><summary>click here for more</summary><p>
https://github.com/raquellondon/Sample_Entropy/blob/main/test_sampen.m
</details>

</details>

## Toolboxes
<details><summary>click here for more</summary><p> 
    
### MATLAB TOOLBOXES 
https://www.physionet.org/content/sampen/1.0.0/

### PYTHON TOOLBOXES
https://github.com/raphaelvallat/antropy </br>
https://cschoel.github.io/nolds/ </br>
https://github.com/MattWillFlood/EntropyHub/tree/main (also MATLAB) </br>

</details>

## Testing SampEn on Raquel's data
<details><summary>click here for more</summary><p>
    
Data from 40 participants, around 800 trials each, 1000 ms pre-stimulus, sr = 1024. </br>

Observations: </br>

Below you can see figures comparing SampEn with the slope of the power spectrum. At first glance we can see the topographies are quite similar, hinting at the idea that these two measures are related. Both have a varying topography and both also vary with time-on-task. What I find weird is that where in the general topography they show an opposite relation (increasing SampEn from posterior to anterior approximately, decreasing PLE), while in the time-on task they go the same way (SampEn decreases, and so does the PLE). This might be something idiosyncratic to my data, but we should keep an eye on this and also think more about the meaning of each measure. 

#### 1.- Average SampEn (averaged over participants and trials) is not evenly distributed across the scalp.
The topography of SampEn (left) is quite similar to the topography of the slope of the power spectrum (right). PLE stands for power law exponent, so that is the x in 1/f^x. A higher x means a steeper slope.

Average SampEn             | Average PLE
:-------------------------:|:-------------------------:
![image](https://github.com/raquellondon/Sample_Entropy/blob/main/Sampen_topo.jpg) | ![image](https://github.com/raquellondon/Sample_Entropy/blob/main/PLEtopo.jpg)

#### 2.- SampEn systematically changes with time on task (remember this is pre-stimulus activity).

Average Correlation SampEn / Trial number             | Average Correlation PLE / Trial number
:-------------------------:|:-------------------------:
![image](https://github.com/raquellondon/Sample_Entropy/blob/main/Time_Sampen_Correlation_topo.jpg) |![image](https://github.com/raquellondon/Sample_Entropy/blob/main/TimeOnTaskSlopeTopo.jpg)
Topography of t-values representing the consistency of the correlation between SampEn and trial order across participants. Electrodes marked in white were significant (FDR correction)|Topography of t-values representing the consistency of the correlation between PLE and trial order across participants. Electrodes marked in white were significant (FDR correction).

Individual Correlation SampEn / Trial number             | Individual Correlation PLE / Trial number
:-------------------------:|:-------------------------:
![image](https://github.com/raquellondon/Sample_Entropy/blob/main/Individual_Correlation_Time_SampEn.jpg) | ![image](https://github.com/raquellondon/Sample_Entropy/blob/main/TimeOnTaskSlopeInd.jpg)
Correlation between SampEn and trial order for each participant averaged across the significant electrodes. | Correlation between PLE and trial order for each participant averaged across the significant electrodes.
</details>

# Open questions:
1.- If r is calculated on each segment separately, then the increase of non-EEG noise (for example) could alter the measurement of SampEn. Should we calulate r on the entire experiment instead of on each trial? I would do it separately for each electrode.</br>
2.- How does the non-stationarity of the data in other aspects influence SampEn? For example, if alpha increases, does that affect our measure directly?</br>




