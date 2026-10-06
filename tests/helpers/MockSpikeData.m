classdef MockSpikeData
% MOCKSPIKEDATA - Stand-in for an AxionFileLoader spike list in allData, for tests
%
% Only GetTimeVoltageVector is provided: Times has one row per sample of each
% waveform cutout (row 1 = cutout start), one column per spike.

    properties
        Times
        Voltages
    end

    methods
        function obj = MockSpikeData(times, voltages)
            obj.Times = times;
            obj.Voltages = voltages;
        end

        function [times, voltages] = GetTimeVoltageVector(obj)
            times = obj.Times;
            voltages = obj.Voltages;
        end
    end
end
