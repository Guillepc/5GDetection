function MonitoreoVivoResourceGridAdaptadoNumCaptures
    %% Configuración inicial y parámetros del receptor SDR y banda 5G
    clear; clc;
    radioOptions = hSDRBase.getDeviceNameOptions;% Obtiene las opciones de dispositivos SDR disponibles en el sistema
    rx = hSDRReceiver(radioOptions{10}); % Crea un objeto receptor SDR usando una de las opciones disponibles, en este caso la opción 10 es B210
    antennaOptions = getAntennaOptions(rx);% Obtiene las opciones de antena para el dispositivo SDR seleccionado
    rx.ChannelMapping = antennaOptions(1); % Establece la antena seleccionada : 1 (RFA:RX2) o 2 (RFB:RX2)
    rx.Gain = 50; % Ajusta la ganancia del receptor SDR
    band = "n78"; % Banda entre 3300-3800 MHz
    GSCN = 7929; % Índice que referencia una frecuencia particular dentro del rango 5G NR para esa banda. 7929 es 5G-laboratorio y 8003 es 5G-calle
    rx.CenterFrequency = hSynchronizationRasterInfo.gscn2frequency(GSCN); % Frecuencia central del SDR usando el GSCN 
    scs = "30kHz"; % Espaciado de subportadora
    nrbSSB = 30;  % Define el número de resource blocks (RBs) ocupados en la SSB 
    scsNumeric = double(extract(scs, digitsPattern)); % % Extrae el valor numérico del espaciado de subportadora
    ofdmInfo = nrOFDMInfo(nrbSSB, scsNumeric); % Calcula información OFDM (incluyendo tasa de muestreo, duración de símbolo, etc)
    rx.SampleRate = ofdmInfo.SampleRate; % Configura la tasa de muestreo del SDR para que coincida con la requerida por OFDM
    
    %% Parámetros de captura y visualización
    interval = 0.3;  % Intervalo entre capturas [s]
    framesPerCapture = 3; % Número de frames por captura. Cada frame corresponde a una duración en tiempo estándar de un frame 5G NR (generalmente 10 ms por frame).
    captureDuration = seconds((framesPerCapture + 1) * 10e-3); % Tiempo total de la duración de cada captura calculado en segundos.Se multiplica framesPerCapture + 1 para asegurarse de que queda un poco más tiempo para capturar toda la ráfaga.
    pauseInterval = 0.00;  % Pausa entre iteraciones
    monitoreoTiempo = 1001 * 100;  % Original tiempo total de monitoreo en segundos aunque con esto se calcula el número de capturas y el tiempo va en funcion a esa variable 
    numCaptures = floor(monitoreoTiempo / interval); % Número de capturas 

    %% Inicializar la figura para mostrar el resource grid
    hFig = figure('Name', 'Monitoreo Vivo Resource Grid 5G Adaptado NumCaptures', 'NumberTitle', 'off');
    hImg = imagesc(zeros(360, 56)); 
    axis xy; colormap('jet'); colorbar;
    xlabel('Símbolos OFDM'); ylabel('Subcarriers');
    title('Resource Grid en vivo 5G');
    hold on; hold off;

    %% Bucle principal de captura y visualización en vivo controlado por numCaptures
    fprintf('Realizando %d capturas con intervalo %.3f s cada una...\n', numCaptures, interval);
    tic;

    for k = 1:numCaptures
        tStart = toc;
        try
            waveform = capture(rx, captureDuration);  % Captura con duration definida
            % Procesar la señal con función findSSBrobusto
            [detectedSSB, gridSSB, burstRect, burstText, ncellid] = findSSBrobusto(waveform, ...
                rx.CenterFrequency, scs, rx.SampleRate);

            % Actualizar imagen solo si detectado y grid válido
            if detectedSSB && ~isempty(gridSSB)
                set(hImg, 'CData', abs(gridSSB));
                title(sprintf('Resource Grid t=%.2fs | CellID=%d | Captura %d/%d', toc, ncellid, k, numCaptures));
                drawnow;
            else
                title(sprintf('No se detectó SSB en t=%.2fs | Captura %d/%d', toc, k, numCaptures));
                drawnow;
            end
        catch ME
            warning('Error en procesamiento en captura %d: %s', k, ME.message);
        end
        
        % Controlar intervalo haciendo pausas para completar el tiempo si es necesario
        tElapsed = toc - tStart;
        tPause = interval - tElapsed;
        if tPause > 0
            pause(tPause);
        end
        % También pausa mínima opcional
        if pauseInterval > 0
            pause(pauseInterval);
        end
    end

    %% Finalizar monitoreo y liberar recursos
    if ishandle(hFig), close(hFig); end
    release(rx);
    fprintf('Monitoreo finalizado tras %d capturas.\n', numCaptures);
end


% Función robusta para encontrar y procesar SSB en la señal capturada

function [detectedSSB, gridSSB, rectSSB, burstText, ncellid] = findSSBrobusto(waveform, centerFrequency, scs, sampleRate)
    % Parámetros de ayuda
    try
        ssbBlockPattern = hSynchronizationRasterInfo.getBlockPattern(scs, centerFrequency);
        scsNumeric = double(extract(scs, digitsPattern));
        searchBW = 3*scsNumeric;
        displayFigure = false;
        
        % Corrección de frecuencia
        [correctedWaveform,~,NID2] = hSSBurstFrequencyCorrect(waveform, ssbBlockPattern, sampleRate, searchBW, displayFigure);
        
        nrbSSB = 20;
        refGrid = zeros([nrbSSB*12 2]);
        refGrid(nrPSSIndices,2) = nrPSS(NID2);
        nSlot = 0;
        
        % Estimación de timing
        timingOffset = nrTimingEstimate(correctedWaveform,nrbSSB,scsNumeric,nSlot,refGrid,SampleRate=sampleRate);
        correctedWaveform = correctedWaveform(1+timingOffset:end,:);
        
        % Demodulación OFDM
        rxGrid = nrOFDMDemodulate(correctedWaveform,nrbSSB,scsNumeric,nSlot,SampleRate=sampleRate);
        
        % Validar tamaño para acceder a símbolo 2:5
        numSymbols = size(rxGrid, 2);
        if numSymbols < 5
            detectedSSB = false;
            gridSSB = [];
            rectSSB = [];
            burstText = '';
            ncellid = -1;
            return;
        end
        
        rxGridSSB = rxGrid(:, 2:5, :); % Símbolos OFDM centrales
        
        % Correlación estimación NID1 usando SSS
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
        
        % Estimar burst SSB más fuerte usando DM-RS
        dmrsIndices = nrPBCHDMRSIndices(ncellid);
        dmrsEst = zeros(1,8);
        for ibar_SSB = 0:7
            refGridDMRS = zeros([240 4]);
            refGridDMRS(dmrsIndices) = nrPBCHDMRS(ncellid,ibar_SSB);
            [hest,nest] = nrChannelEstimate(rxGridSSB,refGridDMRS,'AveragingWindow',[0 1]);
            dmrsEst(ibar_SSB+1) = 10*log10(mean(abs(hest(:).^2)) / nest);
        end
        
        [~, ibar_SSB] = max(dmrsEst);
        strongestSSBIdx = ibar_SSB - 1;
        detectedSSB = maxCorr > 1e-3;
        burstText = sprintf('Strongest SSB: %d', strongestSSBIdx);
        
        % Construir el resource grid para despliegue
        demodRB = 30;
        gridSSB = nrOFDMDemodulate(correctedWaveform, demodRB, scsNumeric, nSlot, SampleRate=sampleRate);
        last = min(56, size(gridSSB, 2));
        gridSSB = gridSSB(:, 1:last, 1);
        
        ssbFreqOrigin = 12*(demodRB-nrbSSB)/2 + 1;
        startSymbol = 1;
        numSymbolsSSB = 4;
        rectSSB = [startSymbol+0.5, ssbFreqOrigin-0.5, numSymbolsSSB, 12*nrbSSB];
        
    catch ME_inner
        % Fallo procesamiento: devolver vacíos y marcado como no detectado
        detectedSSB = false;
        gridSSB = [];
        rectSSB = [];
        burstText = sprintf('Error: %s', ME_inner.message);
        ncellid = -1;
    end
end
