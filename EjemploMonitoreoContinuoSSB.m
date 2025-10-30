% Monitoreo continuo de ráfaga SSB en GSCN 7880 con hSDRReceiver (Wireless Testbench)
clear; clc;

% --- Seleccionar dispositivo SDR ---
radioOptions = hSDRBase.getDeviceNameOptions;
rx = hSDRReceiver(radioOptions{10});           % Ajusta el índice según corresponda
antennaOptions = getAntennaOptions(rx);
rx.ChannelMapping = antennaOptions{1};
rx.Gain = 50;

% --- Configuración del GSCN para banda n77 ---
GSCN = 7880;
freqCenter = hSynchronizationRasterInfo.gscn2frequency(GSCN);
rx.CenterFrequency = freqCenter;

% --- Subcarrier spacing y tasa de muestreo ---
scsOptions = hSynchronizationRasterInfo.getSCSOptions(rx.CenterFrequency);
scs = scsOptions(1);
nrbSSB = 20;
scsNumeric = double(extract(scs, digitsPattern));
ofdmInfo = nrOFDMInfo(nrbSSB, scsNumeric);
rx.SampleRate = ofdmInfo.SampleRate;

% --- Configuración de duración y cantidad de capturas ---
monitorTime = 3;           % Segundos totales
interval = 0.1;             % Intervalo entre capturas (s)
framesPerCapture = 2;
captureDuration = seconds((framesPerCapture+1)*10e-3);
numCaptures = floor(monitorTime / interval);

% --- Variables para métricas ---
timeVec = zeros(numCaptures,1);
powerVec = zeros(numCaptures,1);
snrVec = zeros(numCaptures,1);
detectedVec = false(numCaptures,1);
cellIDVec = zeros(numCaptures,1);
beamIdxVec = zeros(numCaptures,1);
timingOffsetVec = zeros(numCaptures,1);
freqOffsetVec = zeros(numCaptures,1);

disp('Iniciando monitoreo 5G SSB en GSCN 7880...');
tic;

for k = 1:numCaptures
    waveform = capture(rx, captureDuration);

    % Calculo potencia RMS [dB]
    powerVec(k) = 10*log10(mean(abs(waveform).^2));
    % Estimar SNR
    snrVec(k) = estimateSNR(waveform);

    % Detectar SSB y obtener parámetros
    [detectedSSB, ncellid, beamIndex, timingOffset, freqOffset] = ...
        analyzeSSB(waveform, rx.CenterFrequency, scs, rx.SampleRate);

    detectedVec(k) = detectedSSB;
    timeVec(k) = toc;
    if detectedSSB
        cellIDVec(k) = ncellid;
        beamIdxVec(k) = beamIndex;
        timingOffsetVec(k) = timingOffset;
        freqOffsetVec(k) = freqOffset;
        fprintf('[%.1fs] SSB detectado | Pot=%.1f dB | SNR=%.1f dB | CellID=%d | Beam=%d | tOffset=%d | fOffset=%.1f Hz\n', ...
            timeVec(k), powerVec(k), snrVec(k), ncellid, beamIndex, timingOffset, freqOffset);
    else
        fprintf('[%.1fs] No detectado | Pot=%.1f dB | SNR=%.1f dB\n', timeVec(k), powerVec(k), snrVec(k));
    end

    pause(interval);
end

release(rx);

% --- Graficar resultados ---
figure;
subplot(3,1,1);
plot(timeVec, powerVec, 'b.-'); grid on;
xlabel('Tiempo (s)'); ylabel('Potencia (dB)');
title('Evolución de Potencia SSB (GSCN 7880)');

subplot(3,1,2);
plot(timeVec, snrVec, 'r.-'); grid on;
xlabel('Tiempo (s)'); ylabel('SNR (dB)');
title('Evolución de SNR SSB (GSCN 7880)');

subplot(3,1,3);
plot(timeVec(detectedVec), cellIDVec(detectedVec), 'ko-', ...
     timeVec(detectedVec), beamIdxVec(detectedVec), 'm+-');
grid on;
xlabel('Tiempo (s)');
ylabel('CellID / Beam Index');
legend('Cell ID', 'Beam Index');
title('Evolución Cell ID y Beam Index');

% --- Funciones auxiliares ---

function SNRdB = estimateSNR(signal)
    noiseEst = median(abs(signal)).^2;
    signalEst = mean(abs(signal).^2);
    SNRdB = 10*log10(signalEst / noiseEst + eps);
end

function [detectedSSB, ncellid, beamIndex, timingOffset, freqOffset] = analyzeSSB(waveform, centerFrequency, scs, sampleRate)
    % Esta función extiende findSSB para devolver parameters de interés
    
    ssbBlockPattern = hSynchronizationRasterInfo.getBlockPattern(scs,centerFrequency);
    scsNumeric = double(extract(scs,digitsPattern));
    searchBW = 3 * scsNumeric;
    displayFigure = false;

    [correctedWaveform, freqOffset, NID2] = hSSBurstFrequencyCorrect(waveform, ssbBlockPattern, sampleRate, searchBW, displayFigure);

    nrbSSB = 20;
    refGrid = zeros([nrbSSB*12, 2]);
    refGrid(nrPSSIndices, 2) = nrPSS(NID2);

    nSlot = 0;
    timingOffset = nrTimingEstimate(correctedWaveform, nrbSSB, scsNumeric, nSlot, refGrid, SampleRate=sampleRate);
    correctedWaveform = correctedWaveform(1+timingOffset:end, :);

    rxGrid = nrOFDMDemodulate(correctedWaveform, nrbSSB, scsNumeric, nSlot, SampleRate=sampleRate);
    rxGrid = rxGrid(:,2:5,:);

    sssIndices = nrSSSIndices;
    sssRx = nrExtractResources(sssIndices, rxGrid);

    sssEst = zeros(1,336);
    for NID1 = 0:335
        ncellidTemp = (3*NID1) + NID2;
        sssRef = nrSSS(ncellidTemp);
        sssEst(NID1+1) = sum(abs(mean(sssRx .* conj(sssRef),1)).^2);
    end

    [~, idx] = max(sssEst);
    NID1 = idx - 1;
    ncellid = (3*NID1) + NID2;

    dmrsIndices = nrPBCHDMRSIndices(ncellid);
    dmrsEst = zeros(1,8);
    for ibar_SSB = 0:7
        refGridDMRS = zeros([240 4]);
        refGridDMRS(dmrsIndices) = nrPBCHDMRS(ncellid,ibar_SSB);
        [hest,nest] = nrChannelEstimate(rxGrid, refGridDMRS, 'AveragingWindow', [0 1]);
        dmrsEst(ibar_SSB+1) = 10*log10(mean(abs(hest(:).^2)) / nest);
    end

    [~, ibar_SSB] = max(dmrsEst);
    beamIndex = ibar_SSB - 1;

    % Para fines de demostración consideramos detectado si la máxima correlación
    % excede un umbral arbitrario (por ejemplo 1e-3) — este umbral puede ajustarse
    detectedSSB = max(sssEst) > 1e-3;
end
