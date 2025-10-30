clear; clc;

% --- Seleccionar dispositivo SDR ---
radioOptions = hSDRBase.getDeviceNameOptions;
rx = hSDRReceiver(radioOptions{10});           % Ajusta el índice según corresponda
antennaOptions = getAntennaOptions(rx);
rx.ChannelMapping = antennaOptions{1};
rx.Gain = 50;

% --- Configuración del GSCN para banda n77 ---
GSCN = 7929;
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
monitorTime = 8;         % Segundos totales
interval = 0.2;          % Intervalo entre capturas
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
resourceGrids = {};                % Guardar solo los SSB detectados
ssbTimes = [];                     % Vector tiempo para los detectados

disp('Iniciando monitoreo 5G SSB en GSCN 7929...');
tic;

for k = 1:numCaptures
    % Captura de la señal
    waveform = capture(rx, captureDuration);

    % Potencia RMS
    powerVec(k) = 10*log10(mean(abs(waveform).^2));
    % Estimación SNR
    snrVec(k) = estimateSNR(waveform);

    % Análisis SSB ampliado
    [detectedSSB, ncellid, beamIndex, timingOffset, freqOffset, rxGrid] = ...
        analyzeSSB_extended(waveform, rx.CenterFrequency, scs, rx.SampleRate);

    detectedVec(k) = detectedSSB;
    timeVec(k) = toc;

    if detectedSSB
        cellIDVec(k) = ncellid;
        beamIdxVec(k) = beamIndex;
        timingOffsetVec(k) = timingOffset;
        freqOffsetVec(k) = freqOffset;
        resourceGrids{end+1} = abs(rxGrid);  % Solo el valor absoluto
        ssbTimes(end+1) = timeVec(k);
        fprintf('[%.1fs] SSB detectado | Pot=%.1f dB | SNR=%.1f dB | CellID=%d | Beam=%d\n', ...
            timeVec(k), powerVec(k), snrVec(k), ncellid, beamIndex);
    end
    pause(interval);
end

release(rx);

% --- Graficar resultados métricas ---
figure;
subplot(4,1,1);
plot(timeVec, powerVec, 'b.-'); grid on;
xlabel('Tiempo (s)'); ylabel('Potencia (dB)');
title('Evolución de Potencia SSB');

subplot(4,1,2);
plot(timeVec, snrVec, 'r.-'); grid on;
xlabel('Tiempo (s)'); ylabel('SNR (dB)');
title('Evolución SNR SSB');

subplot(4,1,3);
plot(timeVec(detectedVec), cellIDVec(detectedVec), 'ko-', ...
     timeVec(detectedVec), beamIdxVec(detectedVec), 'm+-');
grid on;
xlabel('Tiempo (s)');
ylabel('CellID / Beam');
legend('Cell ID', 'Beam Index');
title('Evolución Cell ID y Beam Index');

% --- Apilar resource grids en el tiempo: heatmap con varias ráfagas ---
subplot(4,1,4);
if ~isempty(resourceGrids)
    % Para simplificar, cortamos todos los resource grids al mismo tamaño (columnas/símbolos)
    minCols = min(cellfun(@(x) size(x,2), resourceGrids)); % Símbolos mínimos entre ráfagas
    minRows = min(cellfun(@(x) size(x,1), resourceGrids)); % Subportadoras mínimas
    % Apilar todos los frames detectados (dimensiones: filas=subportadoras, columnas=símbolos, páginas=tiempo)
    stack = zeros(minRows, minCols, length(resourceGrids));
    for i = 1:length(resourceGrids)
        stack(:,:,i) = resourceGrids{i}(1:minRows,1:minCols);
    end
    % Visualización tipo heatmap: cada fila el "vectorizado" de un resource grid
    waterfallGrid = reshape(stack, minRows*minCols, []); % cada columna = una ráfaga consecutiva
    imagesc(ssbTimes, 1:(minRows*minCols), waterfallGrid);
    xlabel('Tiempo SSB (s)');
    ylabel('Subportadora x Símbolo');
    title('Evolución temporal de Resource Grids (SSB sucesivos)');
    colormap('jet');
    colorbar;
else
    text(0.1,0.5,'No se detectó SSB para graficar resource grids','FontSize',12);
    axis off;
end

%% Funciones auxiliares

function SNRdB = estimateSNR(signal)
    noiseEst = median(abs(signal)).^2;
    signalEst = mean(abs(signal).^2);
    SNRdB = 10*log10(signalEst / noiseEst + eps);
end

function [detectedSSB, ncellid, beamIndex, timingOffset, freqOffset, rxGrid] = analyzeSSB_extended(waveform, centerFrequency, scs, sampleRate)
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

    rxGridSSB = rxGrid(:,2:5,:);
    sssIndices = nrSSSIndices;
    sssRx = nrExtractResources(sssIndices, rxGridSSB);

    sssEst = zeros(1,336);
    for NID1 = 0:335
        ncellidTemp = (3*NID1) + NID2;
        sssRef = nrSSS(ncellidTemp);
        sssEst(NID1+1) = sum(abs(mean(sssRx .* conj(sssRef),1)).^2);
    end

    [maxCorr, idx] = max(sssEst);
    NID1 = idx - 1;
    ncellid = (3*NID1) + NID2;

    dmrsIndices = nrPBCHDMRSIndices(ncellid);
    dmrsEst = zeros(1,8);

    for ibar_SSB = 0:7
        refGridDMRS = zeros([240 4]);
        refGridDMRS(dmrsIndices) = nrPBCHDMRS(ncellid,ibar_SSB);
        [hest,nest] = nrChannelEstimate(rxGridSSB, refGridDMRS, 'AveragingWindow', [0 1]);
        dmrsEst(ibar_SSB+1) = 10*log10(mean(abs(hest(:).^2)) / nest);
    end

    [~, ibar_SSB] = max(dmrsEst);
    beamIndex = ibar_SSB - 1;
    detectedSSB = maxCorr > 1e-3;
end
