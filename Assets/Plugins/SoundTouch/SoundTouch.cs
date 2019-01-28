//----------------------------------------------
//	CreateTime  : #CREATIONDATE#
//	Author      : Taylor Liang
//	Project     : #PROJECTNAME#
//	Company     : #SMARTDEVELOPERS#
//	Instruction : SoundTouch
//	ChangeLog   : None
//----------------------------------------------

using System;
using UnityEngine;

public class SoundTouch
{
    public RateTransposer RateTransposer = new RateTransposer();
    public TDStretch TDStretch = new TDStretch();

    public int channels;
    public int sampleRate;
    public double tempo = 1;
    public double rate = 1;

    private double virtualTempo = 1;
    private double virtualRate = 1;
    private double virtualPitch = 1;

    public void SetSampleRate(int p_sampleRate)
    {
        sampleRate = p_sampleRate;
        TDStretch.SampleRate = p_sampleRate;
    }

    public void SetChannels(int p_channels)
    {
        channels = p_channels;
        RateTransposer.Channels = p_channels;
        TDStretch.Channels = p_channels;
    }

    public void SetTempo(double p_tempo)
    {
        virtualTempo = p_tempo;
        CalcRateAndTempo();
    }

    public void SetRate(double p_rate)
    {
        virtualRate = p_rate;
        CalcRateAndTempo();
    }

    public void SetPitch(double p_pitch)
    {
        virtualPitch = p_pitch;
        CalcRateAndTempo();
    }

    public void SetPitchOctaves(double newPitch)
    {
        virtualPitch = Math.Exp(0.69314718056 * newPitch);
        CalcRateAndTempo();
    }

    public void SetPitchSemiTones(double newPitch)
    {
        SetPitchOctaves(newPitch / 12.0);
    }

    public void CalcRateAndTempo()
    {
        double oldTempo = tempo;
        double oldRate = rate;

        tempo = virtualTempo / virtualPitch;
        rate = virtualPitch * virtualRate;

        RateTransposer.SetRate(rate);
        TDStretch.Tempo = tempo;
    }

    public float[] Process(float[] p_buffer, ref int outputSampleNum)
    {
        TDStretch.SetOverlapLength(8);
        TDStretch.CalcTempoRelated(TDStretch.Tempo);
        TDStretch.CalcSeqParameters();

        float[] result;
        float volOriginal = CalcAvgAmp(p_buffer);
        //float freq = CalcEfficientFrequency(p_buffer, sampleRate, channels);
        //virtualPitch = virtualPitch / (freq / 500);

        if (rate <= 1)
        {
            outputSampleNum = 0;
            result = RateTransposer.Process(p_buffer, ref outputSampleNum);

            SetSampleRate(outputSampleNum);
            result = TDStretch.Process(result);
        }
        else
        {
            result = TDStretch.Process(p_buffer);

            outputSampleNum = 0;
            result = RateTransposer.Process(result, ref outputSampleNum);
        }

        float volResult = CalcAvgAmp(result);
        for (int i = 0; i < result.Length; i++)
        {
            result[i] *= (volOriginal / volResult);
        }

        return result;
    }

    private int avgVolumeWinLength = 1;
    private int maxVolumeWinLength = 1;
    private int frequencyWinLength = 10;

    float CalcEfficientFrequency(float[] p_buffer, int p_sampleRate, int p_channels)
    {
        int winLength = frequencyWinLength;
        float efficientFreq = 0;
        int efficientCount = 0;
        int lastCalcIndex = 0;
        float efficientAmp = CalcMaxAmp(p_buffer) / 3;

        var calcIndex = lastCalcIndex;
        var IsHalfCycle = new Func<int, bool>((index) =>
        {
            if (Math.Abs(p_buffer[calcIndex]) > efficientAmp)
            {
                return (p_buffer[calcIndex] > 0 && p_buffer[index] < 0) ||
                       (p_buffer[calcIndex] < 0 && p_buffer[index] > 0);
            }
            else
            {
                return false;
            }
        });

        for (int i = winLength; i < p_buffer.Length; i += winLength)
        {
            if (Math.Abs(p_buffer[i]) > efficientAmp)
            {
                if (IsHalfCycle(i))
                {
                    float frequency = (float)(i - lastCalcIndex) / p_sampleRate / channels;
                    frequency = 1 / frequency / 2;

                    if (frequency > 200f && frequency < 5000f)
                    {
                        efficientFreq += frequency;
                        efficientCount++;
                    }
                }

                lastCalcIndex = i;
            }
        }

        efficientFreq /= efficientCount;

        return efficientFreq;
    }

    float CalcMaxAmp(float[] p_buffer)
    {
        int winLength = maxVolumeWinLength;
        float volume = 0;

        for (int i = 0; i < p_buffer.Length; i += winLength)
        {
            volume = Math.Max(volume, Math.Abs(p_buffer[i]));
        }

        return volume;
    }

    float CalcAvgAmp(float[] p_buffer)
    {
        int winLength = avgVolumeWinLength;
        float volume = 0;
        int count = 0;

        for (int i = 0; i < p_buffer.Length; i += winLength)
        {
            volume += Math.Abs(p_buffer[i]);
            count++;
        }

        return volume / count;
    }
}