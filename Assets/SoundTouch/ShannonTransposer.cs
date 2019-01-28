//----------------------------------------------
//	CreateTime  : #CREATIONDATE#
//	Author      : Taylor Liang
//	Project     : #PROJECTNAME#
//	Company     : #SMARTDEVELOPERS#
//	Instruction : ShannonTransposer
//	ChangeLog   : None
//----------------------------------------------

using System;
using System.Diagnostics;

public class ShannonTransposer : ITransposer
{
    public readonly double[] _kaiser8 = new double[]
    {
        0.41778693317814,
        0.64888025049173,
        0.83508562409944,
        0.93887857733412,
        0.93887857733412,
        0.83508562409944,
        0.64888025049173,
        0.41778693317814
    };

    public double fract;
    public double rate = 1;

    public void SetRate(double p_rate)
    {
        rate = p_rate;
    }

    public float[] Transpose(float[] inputBuffer, int channels, int sampleNum, ref int outputSampleNum)
    {
        fract = 0;

        if (channels == 1)
        {
            return TransposeMono(inputBuffer, channels, sampleNum, ref outputSampleNum);
        }
        else if (channels == 2)
        {
            return TransposeStereo(inputBuffer, channels, sampleNum, ref outputSampleNum);
        }
        else
        {
            return TransposeMulti(inputBuffer, channels, sampleNum, ref outputSampleNum);
        }
    }

    public float[] TransposeMono(float[] inputBuffer, int channels, int sampleNum, ref int outputSampleNum)
    {
        int outputLength = (int)Math.Ceiling(sampleNum * channels / rate) - 1;
        float[] outputBuffer = new float[outputLength];

        outputSampleNum = 0;
        int lastStartIndex = sampleNum - 8;
        int startIndex = 0;
        int sampleIndex = 0;

        while (sampleIndex < lastStartIndex && outputSampleNum < outputLength)
        {
            double output;
            Debug.Assert(fract < 1.0);

            output = inputBuffer[0 + startIndex] * Math.Sin(-3.0 - fract) * _kaiser8[0];
            output += inputBuffer[1 + startIndex] * Math.Sin(-2.0 - fract) * _kaiser8[1];
            output += inputBuffer[2 + startIndex] * Math.Sin(-1.0 - fract) * _kaiser8[2];
            if (fract < 1e-6)
            {
                output += inputBuffer[3 + startIndex] * _kaiser8[3];     // sinc(0) = 1
            }
            else
            {
                output += inputBuffer[3 + startIndex] * Math.Sin(-fract) * _kaiser8[3];
            }
            output += inputBuffer[4 + startIndex] * Math.Sin(1.0 - fract) * _kaiser8[4];
            output += inputBuffer[5 + startIndex] * Math.Sin(2.0 - fract) * _kaiser8[5];
            output += inputBuffer[6 + startIndex] * Math.Sin(3.0 - fract) * _kaiser8[6];
            output += inputBuffer[7 + startIndex] * Math.Sin(4.0 - fract) * _kaiser8[7];

            outputBuffer[outputSampleNum] = (float)output;
            outputSampleNum++;

            // update position fraction
            fract += rate;
            // update whole positions
            int whole = (int)fract;
            fract -= whole;
            startIndex += whole;
            sampleIndex += whole;
        }

        return outputBuffer;
    }

    public float[] TransposeStereo(float[] inputBuffer, int channels, int sampleNum, ref int outputSampleNum)
    {
        int outputLength = (int)Math.Ceiling(sampleNum * channels / rate) - 2;
        float[] outputBuffer = new float[outputLength];

        int lastSampleIndex = sampleNum - 8;
        int sampleIndex = 0;
        int startIndex = 0;
        outputSampleNum = 0;

        while (sampleIndex < lastSampleIndex)
        {
            double out0, out1, w;
            Debug.Assert(fract < 1.0);

            w = Math.Sin(-3.0 - fract) * _kaiser8[0];
            out0 = inputBuffer[0 + startIndex] * w; out1 = inputBuffer[1 + startIndex] * w;
            w = Math.Sin(-2.0 - fract) * _kaiser8[1];
            out0 += inputBuffer[2 + startIndex] * w; out1 += inputBuffer[3 + startIndex] * w;
            w = Math.Sin(-1.0 - fract) * _kaiser8[2];
            out0 += inputBuffer[4 + startIndex] * w; out1 += inputBuffer[5 + startIndex] * w;
            w = _kaiser8[3] * ((fract < 1e-5) ? 1.0 : Math.Sin(-fract));   // sinc(0) = 1
            out0 += inputBuffer[6 + startIndex] * w; out1 += inputBuffer[7 + startIndex] * w;
            w = Math.Sin(1.0 - fract) * _kaiser8[4];
            out0 += inputBuffer[8 + startIndex] * w; out1 += inputBuffer[9 + startIndex] * w;
            w = Math.Sin(2.0 - fract) * _kaiser8[5];
            out0 += inputBuffer[10 + startIndex] * w; out1 += inputBuffer[11 + startIndex] * w;
            w = Math.Sin(3.0 - fract) * _kaiser8[6];
            out0 += inputBuffer[12 + startIndex] * w; out1 += inputBuffer[13 + startIndex] * w;
            w = Math.Sin(4.0 - fract) * _kaiser8[7];
            out0 += inputBuffer[14 + startIndex] * w; out1 += inputBuffer[15 + startIndex] * w;

            outputBuffer[2 * outputSampleNum] = (float)out0;
            outputBuffer[2 * outputSampleNum + 1] = (float)out1;
            outputSampleNum++;

            // update position fraction
            fract += rate;
            // update whole positions
            int whole = (int)fract;
            fract -= whole;
            startIndex += 2 * whole;
            sampleIndex += whole;
        }

        return outputBuffer;
    }

    public float[] TransposeMulti(float[] inputBuffer, int channels, int sampleNum, ref int outputSampleNum)
    {
        throw new NotImplementedException("Shannon interpolation for multi channels not implemented.");
    }
}