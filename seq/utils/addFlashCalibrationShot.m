function addFlashCalibrationShot(seq, cal, row, rfPhaseRad, adcEvent, labels)
%ADDFLASHCALIBRATIONSHOT Append one RF-spoiled slab FLASH shot of the calibration.
%
%   addFlashCalibrationShot(seq, cal, row, rfPhaseRad, adcEvent, labels)
%
% Inputs:
%   seq         Pulseq sequence object.
%   cal         Calibration designed by designFlashCalibration.
%   row         One acquisition-table row: Wave mode and ky/kz moments [1/m].
%   rfPhaseRad  RF-spoiling phase of this shot [rad].
%   adcEvent    ADC event with its phase, or [] for a dummy shot.
%   labels      Cell array of label events of the ADC ({} for a dummy).
%
% Side effects:
%   Adds the slab excitation, the prephasers, the readout (with the y-sine
%   or z-cosine Wave and the cosine bridges of the row's mode) and the
%   rewinders/spoilers to seq. The shot duration is the calibration TR.

ax = cal.ax;
sys = cal.sys_lowPNS;
rfThis = cal.rf;
rfThis.freqOffset = cal.gzRf.amplitude*cal.slabCenter;
rfThis.phaseOffset = rfPhaseRad-2*pi*rfThis.freqOffset*mr.calcRfCenter(rfThis);

gyPre = mr.makeTrapezoid(ax.d2, 'Area', row.kyArea, ...
    'Duration', cal.preGradDuration, 'system', sys);
gzPre = mr.makeTrapezoid(ax.d3, 'Area', cal.gzRephArea+row.kzArea, ...
    'Duration', cal.preGradDuration, 'system', sys);
gyPost = mr.makeTrapezoid(ax.d2, 'Area', -row.kyArea, ...
    'Duration', cal.spoilDuration, 'system', sys);
gzPost = mr.makeTrapezoid(ax.d3, 'Area', cal.spoilAreaZ-row.kzArea, ...
    'Duration', cal.spoilDuration, 'system', sys);
gxPost = mr.makeTrapezoid(ax.d1, 'Area', cal.groSpoilArea, ...
    'Duration', cal.spoilDuration, 'system', sys);

isWaveY = row.mode == cal.MODE_WAVE_Y;
isWaveZ = row.mode == cal.MODE_WAVE_Z;
readoutEvents = {cal.gxRead};
if isWaveY
    readoutEvents{end+1} = cal.gyWave;
elseif isWaveZ
    readoutEvents{end+1} = cal.gzWave;
end
if ~isempty(adcEvent)
    readoutEvents = [readoutEvents, {adcEvent}, labels];
end

seq.addBlock(rfThis, cal.gzRf, cal.rfBlockEnvelope);
seq.addBlock(mr.align('left', mr.makeDelay(cal.delayTE), ...
    'right', cal.groPre, gyPre, gzPre));
seq.addBlock(cal.gxRampUp);
% The bridge-hold blocks exist only when the calibration carries the z
% Wave; they keep the readout timing of every part identical.
if ~isempty(cal.gxBridgeHold)
    if isWaveZ
        seq.addBlock(cal.gxBridgeHold, cal.gzCosPreBridge);
    else
        seq.addBlock(cal.gxBridgeHold);
    end
end
seq.addBlock(readoutEvents{:});
if ~isempty(cal.gxBridgeHold)
    if isWaveZ
        seq.addBlock(cal.gxBridgeHold, cal.gzCosPostBridge);
    else
        seq.addBlock(cal.gxBridgeHold);
    end
end
seq.addBlock(cal.gxRampDown);
seq.addBlock(gxPost, gyPost, gzPost);
end
