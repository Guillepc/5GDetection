%% SCRIPT PRINCIPAL - MonitoreoFuncionesDisco.m
% Guarda este archivo COMPLETO como "MonitoreoFuncionesDisco.m"
clear; clc; close all;

% === FUNCIONES PRINCIPALES ===
function captureWaveformsDisco(rx, monitorTime, interval, framesPerCapture, sampleRate, outputDir)
    if nargin < 6, outputDir = 'capturas_disco'; end
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end
    
    captureDuration = seconds((framesPerCapture+1)*10e-3);
    numCaptures = floor(monitorTime/interval);
    
    fprintf('=== CAPTURA DIRECTA A DISCO ===\n');
    fprintf('Carpeta: %s/\n', outputDir);
    fprintf('Guardando %d archivos timestamp_*.mat\n', numCaptures);
    
    tic;
    for k = 1:numCaptures
        t_start = tic;
        
        % Captura
        waveform = capture(rx, captureDuration);
        
        % Timestamp único para nombre del archivo
        timestamp = datestr(now, 'yyyymmdd_HHMMSS_FFF');
        filename = fullfile(outputDir, sprintf('timestamp_%s.mat', timestamp));
        
        % Guardar waveform directamente como timestamp_*.mat
        save(filename, 'waveform', 'timestamp', 'k', '-v7.3');
        
        fprintf('[Captura %d/%d] %s (%.3fs) - %.1f MB\n', ...
                k, numCaptures, filename, toc(t_start), ...
                numel(waveform)*8/1e6);
        
        if k < numCaptures
            fprintf('Pausa %.3fs...\n', interval);
            pause(interval);
        end
    end
    fprintf('=== CAPTURA COMPLETADA ===\n');
end

function [resourceGrids, demodTimes, gridSSB1All, ssbTimes, powerVec, snrVec, cellIDVec, ssbRects, ssbTexts] = demodulateAll(waveformsAll, centerFrequency, scs, sampleRate)
    numCaptures = numel(waveformsAll);
    resourceGrids = cell(numCaptures,1);
    demodTimes = zeros(numCaptures,1);
    gridSSB1All = cell(numCaptures,1);
    ssbTimes = zeros(numCaptures,1);
    powerVec = zeros(numCaptures,1);
    snrVec = zeros(numCaptures,1);
    cellIDVec = zeros(numCaptures,1);
    ssbRects = cell(numCaptures,1);
    ssbTexts = cell(numCaptures,1);
    
    fprintf('Demodulando %d waveforms...\n', numCaptures);
    for k = 1:numCaptures
        t_start = tic;
        waveform = waveformsAll{k};
        ssbTimes(k) = k * 0.16;
        
        [~, gridSSB, gridSSB1, rectSSB, burstText, ncellid] = findSSB(waveform, centerFrequency, scs, sampleRate);
        
        powerVec(k) = 10*log10(mean(abs(waveform).^2));
        snrVec(k) = estimateSNR(waveform);
        cellIDVec(k) = ncellid;
        
        resourceGrids{k} = gridSSB;
        gridSSB1All{k} = gridSSB1;
        ssbRects{k} = rectSSB;
        ssbTexts{k} = burstText;
        
        demodTimes(k) = toc(t_start);
        fprintf('[Demod %d] %.3fs | Pot=%.1fdB | SNR=%.1fdB | cellID=%d\n', ...
                k, demodTimes(k), powerVec(k), snrVec(k), ncellid);
    end
end

function gridTimes = generateResourceGrids(resourceGrids, ssbTimes, centerFrequency, ssbRects, ssbTexts)
    t_start = tic;
    visualizeResourceGridsOverlay(resourceGrids, ssbTimes, centerFrequency, ssbRects, ssbTexts);
    gridTimes = toc(t_start);
    fprintf('Tiempo generación Resource Grid: %.3fs\n', gridTimes);
end

% === FUNCIONES AUXILIARES ===
function [detectedSSB, gridSSB, gridSSB1, rectSSB, burstText, ncellid] = findSSB(waveform, centerFrequency, scs, sampleRate)
    ssbBlockPattern = hSynchronizationRasterInfo.getBlockPattern(scs, centerFrequency);
    scsNumeric = double(extract(scs, digitsPattern));
    searchBW = 3*scsNumeric;
    displayFigure = false;
    
    [correctedWaveform, ~, NID2] = hSSBurstFrequencyCorrect(waveform, ssbBlockPattern, sampleRate, searchBW, displayFigure);
    nrbSSB = 20;
    refGrid = zeros([nrbSSB*12 2]);
    refGrid(nrPSSIndices, 2) = nrPSS(NID2);
    nSlot = 0;
    
    timingOffset = nrTimingEstimate(correctedWaveform, nrbSSB, scsNumeric, nSlot, refGrid, SampleRate=sampleRate);
    correctedWaveform = correctedWaveform(1+timingOffset:end, :);
    rxGrid = nrOFDMDemodulate(correctedWaveform, nrbSSB, scsNumeric, nSlot, SampleRate=sampleRate);
    
    rxGridSSB = rxGrid(:, 2:5, :);
    sssIndices = nrSSSIndices;
    sssRx = nrExtractResources(sssIndices, rxGridSSB);
    
    sssEst = zeros(1, 336);
    for NID1 = 0:335
        ncellidTemp = (3*NID1) + NID2;
        sssRef = nrSSS(ncellidTemp);
        sssEst(NID1+1) = sum(abs(mean(sssRx .* conj(sssRef), 1)).^2);
    end
    [maxCorr, idx] = max(sssEst);
    NID1 = idx - 1;
    ncellid = (3*NID1) + NID2;
    
    dmrsIndices = nrPBCHDMRSIndices(ncellid);
    dmrsEst = zeros(1, 8);
    for ibar_SSB = 0:7
        refGrid = zeros([240 4]);
        refGrid(dmrsIndices) = nrPBCHDMRS(ncellid, ibar_SSB);
        [hest, nest] = nrChannelEstimate(rxGridSSB, refGrid, 'AveragingWindow', [0 1]);
        dmrsEst(ibar_SSB+1) = 10*log10(mean(abs(hest(:).^2)) / nest);
    end
    [~, ibar_SSB] = max(dmrsEst);
    strongestSSBIdx = ibar_SSB - 1;
    detectedSSB = maxCorr > 1e-3;
    burstText = sprintf('Strongest SSB: %d', strongestSSBIdx);
    
    demodRB = 45;
    gridSSB1 = nrOFDMDemodulate(correctedWaveform, demodRB, scsNumeric, nSlot, SampleRate=sampleRate);
    last = min(56, size(gridSSB1, 2));
    gridSSB = gridSSB1(:, 1:last, 1);
    ssbFreqOrigin = 12*(demodRB-nrbSSB)/2 + 1;
    startSymbol = 1;
    numSymbolsSSB = 4;
    rectSSB = [startSymbol+0.5 ssbFreqOrigin-0.5 numSymbolsSSB 12*nrbSSB];
end

function visualizeResourceGridsOverlay(resourceGrids, ssbTimes, freqCenter, ssbRects, ssbTexts)
    numCaptures = numel(resourceGrids);
    fixedRows = size(resourceGrids{1}, 1);
    fixedCols = size(resourceGrids{1}, 2);
    allGrids = zeros(fixedRows, fixedCols, numCaptures);
    
    for k = 1:numCaptures
        grid = resourceGrids{k};
        if ~isempty(grid)
            allGrids(1:size(grid,1), 1:size(grid,2), k) = abs(grid);
        end
    end
    
    f = figure('Name', 'Evolución Resource Grid con SSB', 'KeyPressFcn', @keyPressCallback);
    hImage = imagesc(allGrids(:,:,1), 'Parent', gca);
    axis xy; colormap('jet'); colorbar;
    xlabel('OFDM symbol'); ylabel('Subcarrier');
    set(gca, 'XLim', [0.5 fixedCols+0.5], 'YLim', [0.5 fixedRows+0.5]);
    title(sprintf('Resource Grid at t = %.2fs (%.2f MHz) - [1/%d]', ssbTimes(1), freqCenter/1e6, numCaptures));
    
    hSlider = uicontrol('Style', 'slider', 'Min', 1, 'Max', numCaptures, 'Value', 1, ...
        'SliderStep', [1/(numCaptures-1), 1/(numCaptures-1)], 'Position', [150 20 300 20]);
    hSlider.Callback = @(src, ~) updateGrid(round(src.Value));
    
    overlayRect = rectangle('Position', ssbRects{1}, 'EdgeColor', 'r', 'LineWidth', 2);
    overlayText = text(ssbRects{1}(1), ssbRects{1}(2)+20, ssbTexts{1}, 'Color', 'w', 'FontSize', 12);
    
    function updateGrid(idx)
        set(hImage, 'CData', allGrids(:,:,idx));
        set(overlayRect, 'Position', ssbRects{idx});
        set(overlayText, 'Position', [ssbRects{idx}(1), ssbRects{idx}(2)+20]);
        set(overlayText, 'String', ssbTexts{idx});
        title(sprintf('Resource Grid at t = %.2fs (%.2f MHz) - [%d/%d]', ssbTimes(idx), freqCenter/1e6, idx, numCaptures));
        set(gca, 'XLim', [0.5 fixedCols+0.5], 'YLim', [0.5 fixedRows+0.5]);
        drawnow;
    end
    
    function keyPressCallback(~, event)
        idx = round(hSlider.Value);
        switch event.Key
            case 'rightarrow'
                idx = min(idx+1, numCaptures);
            case 'leftarrow'
                idx = max(idx-1, 1);
            otherwise
                return;
        end
        hSlider.Value = idx;
        updateGrid(idx);
    end
end

function SNRdB = estimateSNR(signal)
    noiseEst = median(abs(signal)).^2;
    signalEst = mean(abs(signal).^2);
    SNRdB = 10*log10(signalEst/noiseEst+eps);
end

% =====================================================
% EJECUCIÓN PRINCIPAL
% =====================================================
fprintf('=== INICIANDO MONITOREO CONTINUO A DISCO ===\n');

%% Configuración SDR
radioOptions = hSDRBase.getDeviceNameOptions;
rx = hSDRReceiver(radioOptions{10});
antennaOptions = getAntennaOptions(rx);
rx.ChannelMapping = antennaOptions(1);
rx.Gain = 50;

band = "n78"; GSCN = 7929;
rx.CenterFrequency = hSynchronizationRasterInfo.gscn2frequency(GSCN);
scs = "30kHz";
rx.SampleRate = 19500000;

%% Parámetros PRUEBA (20 capturas)
monitorTime = 1.2;  % ~20 capturas
interval = 0.057;
framesPerCapture = 1;
outputDir = 'capturas_disco';

%% PASO 1: Captura DIRECTA A DISCO
captureWaveformsDisco(rx, monitorTime, interval, framesPerCapture, rx.SampleRate, outputDir);

%% PASO 2: Verificar archivos creados
fprintf('\n=== ARCHIVOS CREADOS EN %s/ ===\n', outputDir);
d = dir(fullfile(outputDir, 'timestamp_*.mat'));
fprintf('Archivos timestamp_*.mat: %d\n', numel(d));
for i = 1:min(5, numel(d))
    fprintf('  %s (%.1f KB)\n', d(i).name, d(i).bytes/1024);
end

%% LIMPIEZA
release(rx);
fprintf('\n=== PROCESO COMPLETADO ===\n');
fprintf('Datos guardados en: %s/\n', outputDir);
fprintf('Archivos: timestamp_YYYYMMDD_HHMMSS_FFF.mat\n');
