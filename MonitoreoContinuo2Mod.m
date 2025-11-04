function [waveformsAll, resourceGrids, ssbTimes, powerVec, snrVec, cellIDVec] = MonitoreoContinuo2Mod
  clear; clc;
  % 1. Configuración SDR, banda y parámetros principales

  radioOptions = hSDRBase.getDeviceNameOptions; % Opciones de dispositivos SDR disponibles en el sistema
  rx = hSDRReceiver(radioOptions{10}); % Selecciona dispositivo SDR (ejemplo: B210)
  antennaOptions = getAntennaOptions(rx); % Opciones de antena para dispositivo seleccionado
  rx.ChannelMapping = antennaOptions(1); % Selecciona antena específica
  rx.Gain = 50; % Ganancia del receptor SDR
 
  band = "n78"; % Banda 5G NR entre 3300-3800 MHz
  GSCN = 8003; % Índice estandarizado de frecuencia para banda seleccionada
  rx.CenterFrequency = hSynchronizationRasterInfo.gscn2frequency(GSCN); % Asigna frecuencia central al SDR

  scs = "30kHz"; % Espaciado subportadora (SCS)
  nrbDemod = 35; % Resource Blocks a demodular, limite ampliado 
  scsNumeric = double(extract(scs, digitsPattern)); % Convierte SCS a valor numérico
  ofdmInfo = nrOFDMInfo(nrbDemod, scsNumeric); % Obtiene info OFDM del sistema (FFT, CP, etc)
  rx.SampleRate = ofdmInfo.SampleRate; % Asigna sample rate acorde a los parámetros OFDM

  %% 2. Parámetros de captura y arreglos para resultados

  monitorTime = 1; % Duración total del monitoreo (en segundos)
  interval = 0.1;  % Intervalo entre ráfagas de captura (en segundos)
  framesPerCapture = 2; % Número de frames por cada ráfaga SDR (cada uno ~10ms)
  captureDuration = seconds((framesPerCapture + 1) * 10e-3); % Tiempo por ráfaga de captura
  
  numCaptures = floor(monitorTime / interval); % Número total de capturas en sesión
  resourceGrids = cell(numCaptures, 1); % Array para guardar los resource grids demodulados
  ssbTimes = zeros(numCaptures, 1); % Tiempos en los que se captura cada ráfaga
  powerVec = zeros(numCaptures, 1); % Vector con potencia media por ráfaga
  snrVec = zeros(numCaptures, 1); % Vector con estimación SNR por ráfaga
  cellIDVec = zeros(numCaptures, 1); % IDs de celda detectados por ráfaga
  waveformsAll = cell(numCaptures, 1); % Cell array para las formas de onda capturadas

  fprintf('Capturando %d ráfagas...\n', numCaptures);
  tic; % Comienza cronómetro real de captura y procesamiento
  

  %% 3. Bucle principal de captura y procesamiento

  for k = 1:numCaptures
      waveform = capture(rx, captureDuration); % Captura burst SDR
      waveformsAll{k} = waveform; % Guarda waveform cruda
      ssbTimes(k) = toc; % Guarda timestamp relativo
      % Detección SSB sencilla, sólo para obtener cellID y stats
      [detectedSSB, ncellid] = detectSSB_simple(waveform, rx.CenterFrequency, scs, rx.SampleRate);
      powerVec(k) = 10 * log10(mean(abs(waveform).^2)); % Potencia media (dB)
      snrVec(k) = estimateSNR(waveform); % Estimación SNR por ráfaga
      cellIDVec(k) = ncellid; % Guarda cellID estimado
      
      %---------------------------------------------
      % Demodulación OFDM para obtener resource grid
      % No se realiza ajuste temporal avanzado aquí
      %---------------------------------------------
      nSlot = 0; % Slot a demodular
      fullGrid = nrOFDMDemodulate(waveform, nrbDemod, scsNumeric, nSlot, "SampleRate", rx.SampleRate); % Demodula grid
      resourceGrids{k} = fullGrid; % Guarda grid completo
      
      % Logging y control visual básico
      fprintf('[%.2fs] Potencia=%.1f dB | SNR=%.1f dB | cellID = %d\n', ...
          ssbTimes(k), powerVec(k), snrVec(k), ncellid);
      pause(interval); % Pausa hasta siguiente captura
  end
  release(rx); % Libera dispositivo SDR para futuras sesiones


  %% 4. Visualización interactiva de todos los resource grids
  visualizeResourceGridsOverlay(resourceGrids, ssbTimes, rx.CenterFrequency);
end

%% Funciones auxilaires

%Función básica para detección rápida de SSB y cellID
function [detectedSSB, ncellid] = detectSSB_simple(waveform, centerFrequency, scs, sampleRate)
  ssbBlockPattern = hSynchronizationRasterInfo.getBlockPattern(scs, centerFrequency); % Patrón de block SSB según banda/frecuencia
  scsNumeric = double(extract(scs, digitsPattern));
  searchBW = 3 * scsNumeric; % BW de búsqueda estándar SSB
  displayFigure = false;
  [correctedWaveform, ~, NID2] = hSSBurstFrequencyCorrect(waveform, ssbBlockPattern, sampleRate, searchBW, displayFigure); % Corrección básica de frecuencia
  nrbSSB = 20; % Resource Blocks para SSB
  refGrid = zeros([nrbSSB * 12, 2]);
  refGrid(nrPSSIndices, 2) = nrPSS(NID2); % Referencia PSS (NID2 estimado)
  nSlot = 0;
  timingOffset = nrTimingEstimate(correctedWaveform, nrbSSB, scsNumeric, nSlot, refGrid, 'SampleRate', sampleRate); % Offset temporal estimado único SSB
  correctedWaveform = correctedWaveform(1 + timingOffset:end, :); % Recorte waveform alineada (SSB)
  rxGrid = nrOFDMDemodulate(correctedWaveform, nrbSSB, scsNumeric, nSlot, 'SampleRate', sampleRate); % Demodulación OFDM básica para SSB
  rxGridSSB = rxGrid(:, 2:5, :); % Extrae símbolo relevante central para SSB
  sssIndices = nrSSSIndices;
  sssRx = nrExtractResources(sssIndices, rxGridSSB); % SSS recibido
  sssEst = zeros(1, 336);
  for NID1 = 0:335
      ncellidTemp = (3 * NID1) + NID2;
      sssRef = nrSSS(ncellidTemp);
      sssEst(NID1 + 1) = sum(abs(mean(sssRx .* conj(sssRef), 1)).^2); % correlación SSS
  end
  [maxCorr, idx] = max(sssEst);
  NID1 = idx - 1;
  ncellid = (3 * NID1) + NID2;
  detectedSSB = maxCorr > 1e-3; % Detección robusta simple
end

%---------------------------------------------------------
% Función para visualizar todos los resource grids
%---------------------------------------------------------
function visualizeResourceGridsOverlay(resourceGrids, ssbTimes, freqCenter)
  numCaptures = numel(resourceGrids);
  fixedRows = size(resourceGrids{1}, 1);
  fixedCols = size(resourceGrids{1}, 2);
  allGrids = zeros(fixedRows, fixedCols, numCaptures); % Matriz 3D para recursos
  for k = 1:numCaptures
      grid = resourceGrids{k};
      if isempty(grid)
          continue;
      end
      allGrids(:, :, k) = abs(grid(:, :, 1)); % Extrae magnitud del grid (slot/antena principal)
  end
  % Gráfica interactiva con slider y control de teclas
  fig = figure('Name', 'Evolución Resource Grid', 'KeyPressFcn', @keyPressCallback);
  hImage = imagesc(allGrids(:, :, 1), 'Parent', gca);
  axis xy;
  colormap('jet');
  colorbar;
  xlabel('OFDM symbol');
  ylabel('Subcarrier');
  set(gca, 'XLim', [0.5, fixedCols + 0.5], 'YLim', [0.5, fixedRows + 0.5]);
  title(sprintf('Resource Grid at t = %.2fs (%.2f MHz) - [1/%d]', ssbTimes(1), freqCenter / 1e6, numCaptures));
  slider = uicontrol('Style', 'slider', 'Min', 1, 'Max', numCaptures, 'Value', 1, ...
      'SliderStep', [1 / (numCaptures - 1), 1 / (numCaptures - 1)], 'Position', [150 20 300 20]);
  slider.Callback = @(src, ~) updateGrid(round(src.Value));
  function updateGrid(idx)
      set(hImage, 'CData', allGrids(:, :, idx));
      title(sprintf('Resource Grid at t = %.2fs (%.2f MHz) - [%d/%d]', ssbTimes(idx), freqCenter / 1e6, idx, numCaptures));
      drawnow;
  end
  function keyPressCallback(~, event)
      idx = round(slider.Value);
      switch event.Key
          case 'rightarrow'
              idx = min(idx + 1, numCaptures);
          case 'leftarrow'
              idx = max(idx - 1, 1);
          otherwise
              return;
      end
      slider.Value = idx;
      updateGrid(idx);
  end
end

%---------------------------------------------------------
% Estimación básica SNR por ráfaga capturada
%---------------------------------------------------------
function SNRdB = estimateSNR(signal)
  noiseEst = median(abs(signal)).^2; % Ruido estimado por mediana
  signalEst = mean(abs(signal).^2); % Señal estimada por media
  SNRdB = 10 * log10(signalEst / noiseEst + eps); % Ratio SNR en dB
end
