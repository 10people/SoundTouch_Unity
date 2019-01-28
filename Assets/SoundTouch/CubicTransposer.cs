//----------------------------------------------
//	CreateTime  : #CREATIONDATE#
//	Author      : Taylor Liang
//	Project     : #PROJECTNAME#
//	Company     : #SMARTDEVELOPERS#
//	Instruction : CubicTransposer
//	ChangeLog   : None
//----------------------------------------------

using System;
using System.Diagnostics;

public class CubicTransposer : ITransposer
{
    public double fract;
    public double rate = 1;

    public readonly float[] _coeffs = new float[] { -0.5f, 1.0f, -0.5f, 0.0f, 1.5f, -2.5f, 0.0f, 1.0f, -1.5f, 2.0f, 0.5f, 0.0f, 0.5f, -0.5f, 0.0f, 0.0f };

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
        int lastStartIndex = sampleNum - 4;
        int startIndex = 0;
        int sampleIndex = 0;

        while (sampleIndex < lastStartIndex)
        {
            float output;
            float x3 = 1.0f;
            float x2 = (float)fract;    // x
            float x1 = x2 * x2;           // x^2
            float x0 = x1 * x2;           // x^3
            float y0, y1, y2, y3;

            Debug.Assert(fract < 1.0);

            y0 = _coeffs[0] * x0 + _coeffs[1] * x1 + _coeffs[2] * x2 + _coeffs[3] * x3;
            y1 = _coeffs[4] * x0 + _coeffs[5] * x1 + _coeffs[6] * x2 + _coeffs[7] * x3;
            y2 = _coeffs[8] * x0 + _coeffs[9] * x1 + _coeffs[10] * x2 + _coeffs[11] * x3;
            y3 = _coeffs[12] * x0 + _coeffs[13] * x1 + _coeffs[14] * x2 + _coeffs[15] * x3;

            output = y0 * inputBuffer[0 + startIndex] + y1 * inputBuffer[1 + startIndex] + y2 * inputBuffer[2 + startIndex] + y3 * inputBuffer[3 + startIndex];

            outputBuffer[outputSampleNum] = output;
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

        int lastSampleIndex = sampleNum - 4;
        int sampleIndex = 0;
        int startIndex = 0;
        outputSampleNum = 0;

        while (sampleIndex < lastSampleIndex)
        {
            float x3 = 1.0f;
            float x2 = (float)fract;    // x
            float x1 = x2 * x2;           // x^2
            float x0 = x1 * x2;           // x^3
            float y0, y1, y2, y3;
            float out0, out1;

            Debug.Assert(fract < 1.0);

            y0 = _coeffs[0] * x0 + _coeffs[1] * x1 + _coeffs[2] * x2 + _coeffs[3] * x3;
            y1 = _coeffs[4] * x0 + _coeffs[5] * x1 + _coeffs[6] * x2 + _coeffs[7] * x3;
            y2 = _coeffs[8] * x0 + _coeffs[9] * x1 + _coeffs[10] * x2 + _coeffs[11] * x3;
            y3 = _coeffs[12] * x0 + _coeffs[13] * x1 + _coeffs[14] * x2 + _coeffs[15] * x3;

            out0 = y0 * inputBuffer[0 + startIndex] + y1 * inputBuffer[2 + startIndex] + y2 * inputBuffer[4 + startIndex] + y3 * inputBuffer[6 + startIndex];
            out1 = y0 * inputBuffer[1 + startIndex] + y1 * inputBuffer[3 + startIndex] + y2 * inputBuffer[5 + startIndex] + y3 * inputBuffer[7 + startIndex];

            outputBuffer[2 * outputSampleNum] = out0;
            outputBuffer[2 * outputSampleNum + 1] = out1;
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

        int lastSampleIndex = sampleNum - 4;
        int sampleIndex = 0;
        int startIndex = 0;
        int outputIndex = 0;
        outputSampleNum = 0;

        while (sampleIndex < lastSampleIndex)
        {
            float x3 = 1.0f;
            float x2 = (float)fract;    // x
            float x1 = x2 * x2;           // x^2
            float x0 = x1 * x2;           // x^3
            float y0, y1, y2, y3;

            Debug.Assert(fract < 1.0);

            y0 = _coeffs[0] * x0 + _coeffs[1] * x1 + _coeffs[2] * x2 + _coeffs[3] * x3;
            y1 = _coeffs[4] * x0 + _coeffs[5] * x1 + _coeffs[6] * x2 + _coeffs[7] * x3;
            y2 = _coeffs[8] * x0 + _coeffs[9] * x1 + _coeffs[10] * x2 + _coeffs[11] * x3;
            y3 = _coeffs[12] * x0 + _coeffs[13] * x1 + _coeffs[14] * x2 + _coeffs[15] * x3;

            for (int c = 0; c < channels; c++)
            {
                outputBuffer[outputIndex] = y0 * inputBuffer[c + startIndex] + y1 * inputBuffer[c + channels + startIndex] + y2 * inputBuffer[c + 2 * channels + startIndex] + y3 * inputBuffer[c + 3 * channels + startIndex];
                outputIndex++;
            }
            outputSampleNum++;

            // update position fraction
            fract += rate;
            // update whole positions
            int whole = (int)fract;
            fract -= whole;
            startIndex += channels * whole;
            sampleNum += whole;
        }

        return outputBuffer;
    }
}