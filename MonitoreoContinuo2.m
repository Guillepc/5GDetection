function [waveformsAll, resourceGrids, ssbTimes, powerVec, snrVec, cellIDVec] = MonitoreoContinuo2
  clear; clc; 

  %  Configuración SDR, banda y parámetros 
    radioOptions = hSDRBase.getDeviceNameOptions; % Obtiene las opciones de dispositivos SDR disponibles en el sistema
    rx = hSDRReceiver(radioOptions{10}); % Crea un objeto receptor SDR usando una de las opciones disponibles, en este caso la opción 10 es B210
    antennaOptions = getAntennaOptions(rx); % Obtiene las opciones de antena para el dispositivo SDR seleccionado
    rx.ChannelMapping = antennaOptions(1); % Establece la antena seleccionada
    rx.Gain = 50; % Ajusta la ganancia del receptor SDR
    % Info de bandas de frecuencia : 
    %fr1BandInfo = hSynchronizationRasterInfo.FR1DLOperatingBand ; 
    band = "n78"; % Banda entre 3300-3800 MHz
    GSCN = 7929; % Índice que referencia una frecuencia particular dentro del rango 5G NR para esa banda
    rx.CenterFrequency = hSynchronizationRasterInfo.gscn2frequency(GSCN); % Frecuencia central del SDR usando el GSCN 
    scs = "30kHz"; % Espaciado de subportadora
    nrbSSB = 20; % Define el número de resource blocks (RBs) ocupados en la SSB 
    scsNumeric = double(extract(scs, digitsPattern)); % % Extrae el valor numérico del espaciado de subportadora
    ofdmInfo = nrOFDMInfo(nrbSSB, scsNumeric); % Calcula información OFDM (incluyendo tasa de muestreo, duración de símbolo, etc)
    rx.SampleRate = ofdmInfo.SampleRate; % Configura la tasa de muestreo del SDR para que coincida con la requerida por OFDM
   
    %% Parámetros de captura

    monitorTime = 1; % Tiempo en el que se va a capturar señales
    interval = 0.1;  % Intervalo entre capturas 
    framesPerCapture = 2; % Número de frames por captura. Cada frame corresponde a una duración en tiempo estándar de un frame 5G NR (generalmente 10 ms por frame).
    captureDuration = seconds((framesPerCapture+1)*10e-3); % Tiempo total de la duración de cada captura calculado en segundos.Se multiplica framesPerCapture + 1 para asegurarse de que queda un poco más tiempo para capturar toda la ráfaga.
    numCaptures = floor(monitorTime/interval); % Número de capturas

    %% Captura secuencial

    resourceGrids = cell(numCaptures,1); % Array para almacenar el resource grid
    ssbTimes = zeros(numCaptures,1); % Vector para guardar el tiempo (en segundos, desde el inicio)
    powerVec = zeros(numCaptures,1); % Vector para almacenar la potencia media (en dB) de cada ráfaga
    snrVec = zeros(numCaptures,1); % Vector para almacenar la estimación de SNR (en dB) de cada ráfaga
    cellIDVec = zeros(numCaptures,1); % Vector para almacenar los identificadores de celda estimados (cellID) para cada ráfaga 
    ssbRects = cell(numCaptures,1); % Celdas para guardar coordenadas (rectángulos) que indican la posición del burst SSB detectado en cada resource grid
    ssbTexts = cell(numCaptures,1); % Celdas para guardar el texto descriptivo (ej. "Strongest SSB: N")
    waveformsAll = cell(numCaptures,1); % Cell array para guardar todas las waveforms

    %% Sección 4: Bucle de captura y procesamiento

    fprintf('Capturando %d ráfagas...\n', numCaptures);
    tic; % Inicia un contador de tiempo para medir cuánto tarda la captura completa

    % Loop principal sobre el número de capturas
    for k = 1:numCaptures
        waveform = capture(rx, captureDuration); % Captura la señal del SDR
        waveformsAll{k} = waveform; % Guardar la forma de onda en el cell array
        ssbTimes(k) = toc;  % Tiempo transcurrido desde que empezó el ciclo

        % --- Detección SSB usando findSSB ---

        [detectedSSB, gridSSB, burstRect, burstText, ncellid] = findSSB(waveform, rx.CenterFrequency, scs, rx.SampleRate); % Detecta SSB más fuerte y obtiene el resource grid
        powerVec(k) = 10*log10(mean(abs(waveform).^2));  % Calcula y almacena la potencia media en dB
        snrVec(k) = estimateSNR(waveform); % Estima el SNR
        cellIDVec(k) = ncellid; % Guarda el ID de celda estimado
        resourceGrids{k} = gridSSB;  % Guarda el resource grid en la celda
        ssbRects{k} = burstRect; % Guarda las coordenadas del burst detectado
        ssbTexts{k} = burstText; % Guarda el texto en la visualización
        fprintf('[%.2fs] Potencia=%.1f dB | SNR=%.1f dB | cellID=%d\n', ssbTimes(k), powerVec(k), snrVec(k), ncellid);
        % Mostrar figura individual estilo ejemplo
              % Mostrar figura individual estilo ejemplo (comentado)
        % figure('Name',sprintf('Resource Grid captura %d',k));
        % imagesc(abs(gridSSB)); axis xy; colormap('jet'); colorbar;
        % xlabel('OFDM symbol'); ylabel('Subcarrier');
        % title(sprintf('Resource Grid en t = %.2fs [%d/%d]', ssbTimes(k), k, numCaptures));
        % if detectedSSB
        %     hold on;
        %     rectangle('Position', burstRect, 'EdgeColor', 'r', 'LineWidth', 2);
        %     text(burstRect(1), burstRect(2)+20, burstText, 'Color', 'w', 'FontSize',12);
        %     hold off;
        % end
        % drawnow;
        pause(interval); % Pausa para mantener el intervalo de captura
    end
    release(rx); % Libera el objeto del SDR para que no quede bloqueado por futuras ejecuciones

%% Visualización final interactiva

    visualizeResourceGridsOverlay(resourceGrids, ssbTimes, rx.CenterFrequency, ssbRects, ssbTexts);

end


%% Funciones auxiliares

% La función findSSB detecta la ráfaga SSB más fuerte en una señal recibida (waveform).
% Realiza corrección de frecuencia, estimación de Timing, demodulación OFDM,
% análisis de secuencias SSS y DM-RS para identificar y marcar la ráfaga.
% Devuelve un booleano indicando si se detectó correcto, el grid de recursos,
% el rectángulo de detección, el texto para visualización y el ID de celda estimado.

function [detectedSSB, gridSSB, rectSSB, burstText, ncellid] = findSSB(waveform,centerFrequency,scs,sampleRate)
    ssbBlockPattern = hSynchronizationRasterInfo.getBlockPattern(scs,centerFrequency);
    scsNumeric = double(extract(scs,digitsPattern));
    searchBW = 3*scsNumeric;
    displayFigure = false;
    [correctedWaveform,~,NID2] = hSSBurstFrequencyCorrect(waveform,ssbBlockPattern,sampleRate,searchBW,displayFigure);
    nrbSSB = 20;
    refGrid = zeros([nrbSSB*12 2]);
    refGrid(nrPSSIndices,2) = nrPSS(NID2);
    nSlot = 0;
    timingOffset = nrTimingEstimate(correctedWaveform,nrbSSB,scsNumeric,nSlot,refGrid,SampleRate=sampleRate);
    correctedWaveform = correctedWaveform(1+timingOffset:end,:);
    rxGrid = nrOFDMDemodulate(correctedWaveform,nrbSSB,scsNumeric,nSlot,SampleRate=sampleRate);

    % Extrae los símbolos centrales donde está el SSB (tipicamente 2:5)
    rxGridSSB = rxGrid(:,2:5,:);
    sssIndices = nrSSSIndices;
    sssRx = nrExtractResources(sssIndices,rxGridSSB);

    % Correlación para estimar NID1 usando SSS
    sssEst = zeros(1,336);
    for NID1 = 0:335
        ncellidTemp = (3*NID1) + NID2;
        sssRef = nrSSS(ncellidTemp);
        sssEst(NID1+1) = sum(abs(mean(sssRx .* conj(sssRef),1)).^2);
    end
    [maxCorr, idx] = max(sssEst);
    NID1 = idx - 1;
    ncellid = (3*NID1) + NID2;
    % Busca el burst SSB más fuerte
    dmrsIndices = nrPBCHDMRSIndices(ncellid);
    dmrsEst = zeros(1,8);
    for ibar_SSB = 0:7
        refGrid = zeros([240 4]);
        refGrid(dmrsIndices) = nrPBCHDMRS(ncellid,ibar_SSB);
        [hest,nest] = nrChannelEstimate(rxGridSSB,refGrid,'AveragingWindow',[0 1]);
        dmrsEst(ibar_SSB+1) = 10*log10(mean(abs(hest(:).^2)) / nest);
    end
    [~, ibar_SSB] = max(dmrsEst);
    strongestSSBIdx = ibar_SSB - 1;
    detectedSSB = maxCorr > 1e-3;
    burstText = sprintf('Strongest SSB: %d', strongestSSBIdx);

    % RECURSO GRID Y RECTÁNGULO ESTILO EJEMPLO
    demodRB = 30;
    gridSSB = nrOFDMDemodulate(correctedWaveform,demodRB,scsNumeric,nSlot,SampleRate=sampleRate);
    last = min(56,size(gridSSB,2));
    gridSSB = gridSSB(:,1:last,1);
    ssbFreqOrigin = 12*(demodRB-nrbSSB)/2 + 1;
    startSymbol = 1;
    numSymbolsSSB = 4;
    rectSSB = [startSymbol+0.5 ssbFreqOrigin-0.5 numSymbolsSSB 12*nrbSSB];
end

% La función visualizeResourceGridsOverlay genera una figura interactiva para visualizar la evolución
% del resource grid a lo largo de múltiples capturas. Incluye un slider que
% permite navegar por los distintos instantes y muestra un rectángulo y texto
% que marcan la posición del burst SSB más fuerte correspondiente a cada momento.

function visualizeResourceGridsOverlay(resourceGrids, ssbTimes, freqCenter, ssbRects, ssbTexts)
    numCaptures = numel(resourceGrids);
    fixedRows = size(resourceGrids{1},1);
    fixedCols = size(resourceGrids{1},2);
    allGrids = zeros(fixedRows, fixedCols, numCaptures);
    for k = 1:numCaptures
        grid = resourceGrids{k};
        if isempty(grid), continue; end
        allGrids(1:size(grid,1), 1:size(grid,2), k) = abs(grid);
    end
    f = figure('Name','Evolución Resource Grid con SSB','KeyPressFcn',@keyPressCallback);
    hImage = imagesc(allGrids(:,:,1), 'Parent', gca);
    axis xy; colormap('jet'); colorbar;
    xlabel('OFDM symbol'); ylabel('Subcarrier');
    set(gca, 'XLim', [0.5 fixedCols+0.5], 'YLim', [0.5 fixedRows+0.5]);
    title(sprintf('Resource Grid at t = %.2fs (%.2f MHz) - [1/%d]', ssbTimes(1), freqCenter/1e6, numCaptures));
    hSlider = uicontrol('Style','slider','Min',1,'Max',numCaptures,'Value',1,...
        'SliderStep', [1/(numCaptures-1), 1/(numCaptures-1)],'Position', [150 20 300 20]);
    hSlider.Callback = @(src,~) updateGrid(round(src.Value));
    overlayRect = rectangle('Position', ssbRects{1}, 'EdgeColor', 'r', 'LineWidth', 2);
    overlayText = text(ssbRects{1}(1), ssbRects{1}(2)+20, ssbTexts{1}, 'Color','w','FontSize',12);
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


% La función estimate SNR calcula y devuelve la estimación del Signal-to-Noise Ratio (SNR) en dB
% para la señal recibida, usando la mediana para estimar el ruido y la media
% para la señal, proporcionando una medida de la calidad de la señal recibida.

function SNRdB = estimateSNR(signal)
    noiseEst = median(abs(signal)).^2;
    signalEst = mean(abs(signal).^2);
    SNRdB = 10*log10(signalEst/noiseEst+eps);
end
