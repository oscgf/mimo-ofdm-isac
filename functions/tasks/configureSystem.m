%% Subfunción: configuración del sistema
function [systemParams, transmitter, receiver, antennas, arrays] = configureSystem(params)
    % This function configures the system parameters, transmitter, receiver,
    % antennas, and arrays based on the provided parameters.

    % Initialize structures
    systemParams = struct();
    transmitter = struct();
    receiver = struct();
    antennas = struct();
    arrays = struct();

    % Initial parameters
    systemParams.carrierFrequency = params.CarrierFreq_GHz * 1e9; % Carrier frequency (Hz)
    systemParams.waveLength = freq2wavelen(systemParams.carrierFrequency); % Wavelength (m)
    systemParams.bandwidth = params.Bandwidth_MHz * 1e6; % Bandwidth (Hz)
    systemParams.sampleRate = systemParams.bandwidth; % Sample rate (Hz)
    
    systemParams.NumSubcarriers = params.NumSubcarriers;  % Number of subcarriers
    systemParams.numDataStreams = params.numDataStreams; % Number of data streams
    systemParams.Nframe = params.Nframe;                 % Total number of transmitted OFDM frames
    systemParams.bitsPerSymbol = params.bitsPerSymbol;   % Bits per QAM symbol (and OFDM data subcarrier)
    systemParams.modOrder = 2^systemParams.bitsPerSymbol;        % Modulation order

    % Transmitter
    % - Parameters
    systemParams.peakPower = params.peakPower_W; % Peak power (W)
    % - Transmitter configuration
    transmitter = phased.Transmitter(...
        'PeakPower', systemParams.peakPower, ...
        'Gain', 0);

    % Receiver
    % - Parameters
    systemParams.noiseFigure = params.noiseFigure_dB; % Noise figure (dB)
    systemParams.referenceTemperature = params.referenceTemperature_K; % Reference temperature (K)
    % - Receiver configuration
    receiver = phased.Receiver(...
        'SampleRate', systemParams.sampleRate,...
        'NoiseFigure', systemParams.noiseFigure,...
        'ReferenceTemperature', systemParams.referenceTemperature,...
        'AddInputNoise', true,...
        'InputNoiseTemperature', systemParams.referenceTemperature,...
        'Gain', 0);

    % Antenna Parameters (ULA)
    antennas.Ntx = params.TxAntennas; % Number of transmit antennas
    antennas.Nrx = params.RxAntennas; % Number of receive antennas

    element = phased.IsotropicAntennaElement('BackBaffled', true); % Isotropic antenna element with back baffling

    arrays.tx = phased.ULA(antennas.Ntx, systemParams.waveLength/2, 'Element', element); % Transmitter array
    arrays.rx = phased.ULA(antennas.Nrx, systemParams.waveLength/2, 'Element', element); % Receiver array
end