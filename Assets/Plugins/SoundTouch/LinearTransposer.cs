//----------------------------------------------
//	CreateTime  : #CREATIONDATE#
//	Author      : Taylor Liang
//	Project     : #PROJECTNAME#
//	Company     : #SMARTDEVELOPERS#
//	Instruction : LinearTransposer
//	ChangeLog   : None
//----------------------------------------------

using System;
using System.Diagnostics;

public class LinearTransposer : ITransposer
{
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
        int lastStartIndex = sampleNum - 1;
        int startIndex = 0;
        int sampleIndex = 0;

        while (sampleIndex < lastStartIndex && outputSampleNum < outputLength)
        {
            Debug.Assert(fract < 1.0);

            double output = (1.0 - fract) * inputBuffer[startIndex] + fract * inputBuffer[startIndex + 1];
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

        int lastSampleIndex = sampleNum - 1;
        int sampleIndex = 0;
        int startIndex = 0;
        outputSampleNum = 0;

        while (sampleIndex < lastSampleIndex)
        {
            double out0, out1;
            Debug.Assert(fract < 1.0);

            out0 = (1.0 - fract) * inputBuffer[startIndex] + fract * inputBuffer[startIndex + 2];
            out1 = (1.0 - fract) * inputBuffer[startIndex + 1] + fract * inputBuffer[startIndex + 3];
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
        int outputLength = (int)Math.Ceiling(sampleNum * channels / rate) - channels;
        float[] outputBuffer = new float[outputLength];

        int lastSampleIndex = sampleNum - 1;
        int sampleIndex = 0;
        int startIndex = 0;
        int outputIndex = 0;
        outputSampleNum = 0;

        while (sampleIndex < lastSampleIndex)
        {
            float temp, vol1, fract_float;

            vol1 = (float)(1.0 - fract);
            fract_float = (float)fract;
            for (int c = 0; c < channels; c++)
            {
                temp = vol1 * inputBuffer[startIndex + c] + fract_float * inputBuffer[startIndex + c + channels];
                outputBuffer[outputIndex] = (float)temp;
                outputIndex++;
            }
            outputSampleNum++;

            fract += rate;

            int iWhole = (int)fract;
            fract -= iWhole;
            sampleIndex += iWhole;
            startIndex += iWhole * channels;
        }

        return outputBuffer;
    }
}