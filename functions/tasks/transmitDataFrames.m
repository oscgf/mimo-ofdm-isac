%% Subfunción: transmisión de tramas OFDM y evaluación de BER
function [radarDataCube] = transmitDataFrames(systemParams, transmitter, receiver, antennas, scenario, channel, ofdm, precoding)
    % This function transmits OFDM frames and evaluates the Bit Error Rate (BER) 

    % Subframe A
    [subframeAMod, subframeADemod, subframeAInfo] = createModDemod(true, ofdm, antennas);

    % Subframe B
    % Define pilot subcarrier indices for subframe B
    [ofdm.pilotIdxs, ofdm.pilots] = generatePilots(ofdm, antennas.Ntx);
    [subframeBMod, subframeBDemod, subframeBInfo] = createModDemod(false, ofdm, antennas);

    % Indices of data subcarriers in the subframe B
    ofdm.subframeBdataSubcarrierIdxs = setdiff(ofdm.numGuardBandCarriers(1)+1:(ofdm.Nsub - ofdm.numGuardBandCarriers(2)), ofdm.pilotIdxs);

    fprintf("Velocity resolution: %.2f (m/s).\n", dop2speed(1/(ofdm.Nframe*ofdm.Tofdm*ofdm.Mt), systemParams.waveLength));

    % Input data size for subframe A and B
    subAInputSize = [subframeAInfo.DataInputSize(1), subframeAInfo.DataInputSize(2), ofdm.numDataStreams];
    subBInputSize = [subframeBInfo.DataInputSize(1), subframeBInfo.DataInputSize(2), ofdm.numDataStreams];

    % Initialize radar data cube
    radarDataCube = zeros(ofdm.numActiveSubcarriers, antennas.Nrx, ofdm.Nframe);

    % Simulate formation, transmission, and reception of an OFDM frame one at a time.
    for i = 1:ofdm.Nframe
        % Generate binary payload for subframes A and B and modulate data using QAM
        [subframeABin, subA] = generateSubframe(true, subframeAMod, subAInputSize, systemParams, ofdm, precoding, antennas);
        [subframeBBin, subB] = generateSubframe(false, subframeBMod, subBInputSize, systemParams, ofdm, precoding, antennas);

        % Binary data transmitted in the ith frame
        txDataBin = cat(1, subframeABin(:), subframeBBin(:));

        % Reshape and combine subframes A and B to transmit the whole frame one symbol at a time
        subA = reshape(subA, ofdm.ofdmSymbolLengthWithCP, ofdm.subframeALength, []);
        subB = reshape(subB, ofdm.ofdmSymbolLengthWithCP, antennas.Ntx, []);
        ofdmSignal = [subA subB];

        % Preallocate space for the received signal
        rxSignal = zeros(size(ofdmSignal,1), size(ofdmSignal,2), antennas.Nrx);

        % Transmit one OFDM symbol at a time
        for s = 1:size(ofdmSignal,2)
            % Update target positions
            [scenario.targetPositions, scenario.targetVelocities] = scenario.targetMotion(ofdm.Tofdm);

            % Transmit signal
            tx = transmitter(squeeze(ofdmSignal(:,s,:)));

            % Apply scattering MIMO channel propagation effects
            chanOut = channel(tx,...
                [scenario.scatterPos scenario.targetPositions], ...
                [zeros(size(scenario.scatterPos)) scenario.targetVelocities], ...
                [scenario.scatterRC scenario.targetRC]);

            % Add thermal noise at the receiver
            rxSignal(:,s,:) = receiver(chanOut);
        end

        % Separate the received signal into subframes A and B
        rxSubframeA = rxSignal(:, 1:ofdm.subframeALength, :);
        rxSubframeA = reshape(rxSubframeA, [], antennas.Nrx);

        rxSubframeB = rxSignal(:, ofdm.subframeALength+1:end, :);
        rxSubframeB = reshape(rxSubframeB, [], antennas.Nrx);

        % Demodulate subframe A and B and apply the combining weights
        [rxSubframeAQamComb, ~] = demodulateAndApplyWeights(true, subframeADemod, rxSubframeA, ofdm, precoding, antennas);
        [rxSubframeBQamComb, rxPilots] = demodulateAndApplyWeights(false, subframeBDemod, rxSubframeB, ofdm, precoding, antennas);

        % Demodulate the QAM data and compute the bit error rate for the ith frame
        rxDataQam = cat(1, rxSubframeAQamComb(:), rxSubframeBQamComb(:));
        rxDataBin = qamdemod(rxDataQam, ofdm.modOrder, 'OutputType', 'bit', 'UnitAveragePower', true);
        [~, ratio] = biterr(txDataBin, rxDataBin);
        fprintf("Frame %d bit error rate: %.4f\n", i, ratio);

        % Update channel
        %  - Estimate channel matrix using pilots in the subframe B
        channelMatrix = helperInterpolateChannelMatrix(ofdm.Nsub, ofdm.numGuardBandCarriers, ofdm.pilots, rxPilots, ofdm.pilotIdxs);

        %  - Compute precoding and combining weights for the next frame
        [Wp, Wc, ~, G] = diagbfweights(channelMatrix);
        precoding.Wp = Wp;
        precoding.Wc = Wc;
        precoding.G = G;

        % Store the radar data
        radarDataCube(:,:,i) = squeeze(sum(channelMatrix,2));
    end

    % Constellation diagram setup
    refconst = qammod(0:ofdm.modOrder-1, ofdm.modOrder, 'UnitAveragePower', true);
    constellationDiagram = comm.ConstellationDiagram('NumInputPorts', 1, ...
        'ReferenceConstellation', refconst, 'ChannelNames', {'Received QAM Symbols'});

    constellationDiagram(rxDataQam);
end
