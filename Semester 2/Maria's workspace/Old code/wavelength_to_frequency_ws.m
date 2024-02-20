
function frequencies = wavelength_to_frequency_ws(lambda_min, lambda_max)
    % Constants
    c = 299792458; % Speed of light in meters per second

    % Calculate corresponding frequencies
    f_min = c / lambda_max; % Minimum frequency in Hz
    f_max = c / lambda_min; % Maximum frequency in Hz

    % Generate array of frequencies
    frequencies = linspace(f_min, f_max, 102); % Adjust the number of points as needed
end

