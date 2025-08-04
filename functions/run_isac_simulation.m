function run_isac_simulation(params)
    % Set the random number generator for reproducibility
    rng('default');

    %% 1. System Parameters
    [systemParams, transmitter, receiver, antennas, arrays] = configureSystem(params);

    %% 2. Escenario ISAC
    [scenario, channel] = configureScenario(systemParams, arrays);

    %% 3. Configuración OFDM y estimación de canal (initial channel sounding)
    ofdm = configureOFDM(systemParams, antennas, scenario);
    [~, precoding] = initialChannelEstimation(systemParams, transmitter, receiver, antennas, arrays, scenario, channel, ofdm);

    %% 4. Transmisión de tramas OFDM y evaluación de BER
    profile on
    [radarDataCube] = transmitDataFrames(systemParams, transmitter, receiver, antennas, scenario, channel, ofdm, precoding);
    profile viewer
    
    %% 5. Procesado radar y métricas de sensado.
    processRadarData(radarDataCube, systemParams, scenario, arrays, ofdm);
end