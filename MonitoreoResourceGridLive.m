function MonitoreoVivoResourceGrid
 %% Configuración inicial y parámetros del receptor SDR y banda 5G

    clear; clc; % Limpia el espacio de trabajo y la ventana de comandos para empezar fresco.
    radioOptions = hSDRBase.getDeviceNameOptions; % Obtiene una lista de dispositivos SDR disponibles en el sistema.
    rx = hSDRReceiver(radioOptions{10}); % Crea un objeto receptor SDR usando el dispositivo listado en la posición 10 (B210).
    antennaOptions = getAntennaOptions(rx); % Obtiene las opciones de antena del receptor
    rx.ChannelMapping = antennaOptions(1); % Configura el receptor para usar la primera opción de antena disponible: 1 (RFA-RX2) y 2 (RFB-RX2)
    rx.Gain = 50; % Configura la ganancia del receptor
    band = "n78"; % Banda entre 3300-3800 MHz
    GSCN = 8003; % Índice que referencia  una frecuencia particular dentro del rango 5G NR para esa banda:7929(5G-Laboratorio), 8003/7880(5G-Calle)
    rx.CenterFrequency = hSynchronizationRasterInfo.gscn2frequency(GSCN); % Convierte GSCN a frecuencia central y la asigna al receptor SDR.
    scs = "30kHz"; % Define el subcarrier spacing
    nrbSSB = 30; % Número de Resource Blocks
    scsNumeric = double(extract(scs,digitsPattern)); % Extrae el valor numérico 30 de la cadena "30kHz"
    ofdmInfo = nrOFDMInfo(nrbSSB, scsNumeric); % Obtiene información OFDM
    rx.SampleRate = ofdmInfo.SampleRate; % Configura la tasa de muestreo del receptor con el valor calculado para OFDM.
    
 %% Parámetros de captura y visualización
    interval = 0.145;  % Intervalo entre capturas
    framesPerCapture = 3; % Número de frames por captura. Cada frame corresponde a una duración en tiempo estándar de un frame 5G NR (generalmente 10 ms por frame).
    captureDuration = seconds((framesPerCapture+1)*10e-3); % Tiempo total de la duración de cada captura calculado en segundos.Se multiplica framesPerCapture + 1 para asegurarse de que queda un poco más tiempo para capturar toda la ráfaga.  % Duración de cada captura (~3 frames de 10ms cada uno)
    pauseInterval = 0.00;              % Pausa entre iteraciones en segundos
    monitoreoTiempo = 600*2;          % Duración total monitorización en segundos (e.g. 20 min)
    %numCaptures = floor(monitorTime/interval); % Número de capturas
    %REDEFINIR CAPTURE CON EL NUMCAPTURES NO CON CAPTURE DURATION
    
 %% Inicializar la figura para mostrar el resource grid    
 
    hFig = figure('Name','Monitoreo Vivo Resource Grid 5G','NumberTitle','off'); % Crea una ventana de figura con título personalizado y sin mostrar número de figura.
    hImg = imagesc(zeros(360,56)); % Ajustar tamaño según configuración (filas=subcarriers, cols=símbolos OFDM)
    axis xy; colormap('jet'); colorbar; % Configura el eje para que el origen esté en la esquina inferior izquierda
    xlabel('Símbolos OFDM'); ylabel('Subcarriers'); % Aplica la paleta de colores 'jet'
    title('Resource Grid en vivo 5G');
    hold on; % Prepara la figura para posibles posteriores gráficas sin borrar la actual 
    hold off;
    
 %% Iniciar temporizador
    
    tic; % Inicia un contador para medir el tiempo de ejecución del monitoreo.

    
 %% Bucle principal de captura y visualización en vivo
    
    while ishandle(hFig) && toc < monitoreoTiempo  % El bucle se ejecuta mientras la figura está abierta y el tiempo transcurrido sea menor que el total de monitoreo.

        try
            waveform = capture(rx, captureDuration);  % Captura una señal del receptor SDR durante la duración definida.
            % Procesar la señal con función findSSB adaptada
            [detectedSSB, gridSSB, burstRect, burstText, ncellid] = findSSBrobusto(waveform,...
            rx.CenterFrequency, scs, rx.SampleRate);
            
            % Solo actualizar si grid válido y detectado
            if detectedSSB && ~isempty(gridSSB) % Actualiza la imagen con el resource grid detectado
                set(hImg, 'CData', abs(gridSSB));
                
          
                title(sprintf('Resource Grid t=%.2fs | CellID=%d', toc, ncellid)); % Actualiza el título con tiempo transcurrido y Cell ID detectado
                drawnow;
            else
                % Opcional: mostrar mensaje si no se detecta SSB
                title(sprintf('No se detectó SSB en t=%.2fs', toc));
                drawnow;
            end
        catch ME
            warning('Error en procesamiento: %s', ME.message); % Muestra advertencia si ocurre algún error en esta iteración.
            % Continuar sin detener la captura
        end
        pause(pauseInterval); % Controla la tasa de actualización
    end
    
 %% Finalizar monitoreo y liberar recursos

    if ishandle(hFig), close(hFig); end
    release(rx);
    fprintf('Monitoreo finalizado.\n');
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
