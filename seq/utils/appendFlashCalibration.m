function appendFlashCalibration(seq, cal)
%APPENDFLASHCALIBRATION Append the designed 3D FLASH calibration to a sequence.
%
%   appendFlashCalibration(seq, cal)
%
% Inputs:
%   seq  Pulseq sequence object; the calibration blocks are added at its end.
%   cal  Calibration designed by designFlashCalibration.
%
% Side effects:
%   Adds the dummy and acquired shots of every calibration part to seq. Each
%   part starts with its dummy or settling shots, and the RF-spoiling phase
%   restarts at the first acquired shot, so that a Wave part and its
%   reference share the same RF/ADC phase evolution. Every ADC carries its
%   compact LIN label, zero PAR/SLC/SET/ECO and REF=true.

for p = 1:numel(cal.parts)
    rows = cal.acqTable(cal.partStart(p):cal.partStop(p));
    % Dummies cycle through the last lines of the part so that the encoding
    % sequence entering the first acquired shot is the acquired one.
    dummyRows = mod((-cal.parts(p).nDummy:-1), numel(rows))+1;
    rfPhase = 0;
    rfIncrement = 0;
    for iDummy = 1:numel(dummyRows)
        addFlashCalibrationShot(seq, cal, rows(dummyRows(iDummy)), ...
            rfPhase/180*pi, [], {});
        rfIncrement = mod(rfIncrement+cal.rfSpoilingInc, 360);
        rfPhase = mod(rfPhase+rfIncrement, 360);
    end

    rfPhase = 0;
    rfIncrement = 0;
    for iRow = 1:numel(rows)
        row = rows(iRow);
        adcThis = cal.adc;
        adcThis.phaseOffset = rfPhase/180*pi;
        labels = [{cal.lblLIN{row.rawLin0+1}}, cal.lblZero];
        addFlashCalibrationShot(seq, cal, row, rfPhase/180*pi, adcThis, labels);
        rfIncrement = mod(rfIncrement+cal.rfSpoilingInc, 360);
        rfPhase = mod(rfPhase+rfIncrement, 360);
    end
end
end
