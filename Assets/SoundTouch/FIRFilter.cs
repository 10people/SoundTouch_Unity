//----------------------------------------------
//	CreateTime  : 14/01/2019 16:20:11
//	Author      : Taylor Liang
//	Project     : My_Demo
//	Company     : BuGu
//	Instruction : FIRFilter
//	ChangeLog   : None
//----------------------------------------------

using System;
using UnityEngine;

public class FIRFilter
{
    // Number of FIR filter taps
    int length;
    // Number of FIR filter taps divided by 8
    int lengthDiv8;

    // Result divider factor in 2^k format
    int resultDivFactor;

    // Result divider value.
    float resultDivider;

    // Memory for filter coefficients
    float[] filterCoeffs;

    // Set filter coeffiecients and length.
    //
    // Throws an exception if filter length isn't divisible by 8
    public void SetCoefficients(float[] coeffs, int newLength, int uResultDivFactor)
    {
        Debug.Assert(newLength > 0);

        if (newLength % 8 != 0)
        {
            Debug.LogError("FIR filter length not divisible by 8");
        }

        lengthDiv8 = newLength / 8;
        length = lengthDiv8 * 8;
        Debug.Assert(length == newLength);

        resultDivFactor = uResultDivFactor;
        resultDivider = (float)Math.Pow(2.0, (int)resultDivFactor);

        filterCoeffs = new float[length];
        Array.Copy(coeffs, filterCoeffs, length);
    }

    // Applies the filter to the given sequence of samples. 
    //
    // Note : The amount of outputted samples is by value of 'filter_length' 
    // smaller than the amount of input samples.
    public float[] Evaluate(float[] src, int numSamples, int numChannels, ref int outputSampleNum)
    {
        Debug.Assert(length > 0);
        Debug.Assert(lengthDiv8 * 8 == length);

        if (numSamples < length)
        {
            outputSampleNum = 0;
            return src;
        }

        if (numChannels == 1)
        {
            return EvaluateFilterMono(src, numSamples, ref outputSampleNum);
        }
        else if (numChannels == 2)
        {
            return EvaluateFilterStereo(src, numSamples, ref outputSampleNum);
        }
        else
        {
            Debug.Assert(numChannels > 0);
            return EvaluateFilterMulti(src, numSamples, numChannels, ref outputSampleNum);
        }
    }

    // Usual C-version of the filter routine for stereo sound
    float[] EvaluateFilterStereo(float[] src, int numSamples, ref int outputSampleNum)
    {
        int j, end;
        // when using floating point samples, use a scaler instead of a divider
        // because division is much slower operation than multiplying.
        double dScaler = 1.0 / (double)resultDivider;

        Debug.Assert(length != 0);
        Debug.Assert(src.Length > 0);
        Debug.Assert(filterCoeffs.Length > 0);

        end = 2 * (numSamples - length);
        outputSampleNum = numSamples - length;
        float[] dest = new float[end + 1];

        for (j = 0; j < end; j += 2)
        {
            double suml, sumr;
            int i;

            suml = sumr = 0;

            for (i = 0; i < length; i += 4)
            {
                // loop is unrolled by factor of 4 here for efficiency
                suml += src[j + 2 * i + 0] * filterCoeffs[i + 0] +
                        src[j + 2 * i + 2] * filterCoeffs[i + 1] +
                        src[j + 2 * i + 4] * filterCoeffs[i + 2] +
                        src[j + 2 * i + 6] * filterCoeffs[i + 3];
                sumr += src[j + 2 * i + 1] * filterCoeffs[i + 0] +
                        src[j + 2 * i + 3] * filterCoeffs[i + 1] +
                        src[j + 2 * i + 5] * filterCoeffs[i + 2] +
                        src[j + 2 * i + 7] * filterCoeffs[i + 3];
            }

            suml *= dScaler;
            sumr *= dScaler;
            dest[j] = (float)suml;
            dest[j + 1] = (float)sumr;
        }
        return dest;
    }


    // Usual C-version of the filter routine for mono sound
    float[] EvaluateFilterMono(float[] src, int numSamples, ref int outputSampleNum)
    {
        int j;
        // when using floating point samples, use a scaler instead of a divider
        // because division is much slower operation than multiplying.
        double dScaler = 1.0 / (double)resultDivider;

        Debug.Assert(length != 0);

        outputSampleNum = numSamples - length;
        float[] dest = new float[outputSampleNum];

        for (j = 0; j < outputSampleNum; j++)
        {
            double sum;
            int i;

            sum = 0;
            for (i = 0; i < length; i += 4)
            {
                // loop is unrolled by factor of 4 here for efficiency
                sum += src[j + i + 0] * filterCoeffs[i + 0] +
                       src[j + i + 1] * filterCoeffs[i + 1] +
                       src[j + i + 2] * filterCoeffs[i + 2] +
                       src[j + i + 3] * filterCoeffs[i + 3];
            }
            sum *= dScaler;
            dest[j] = (float)sum;
        }
        return dest;
    }


    float[] EvaluateFilterMulti(float[] src, int numSamples, int numChannels, ref int outputSampleNum)
    {
        int j, end;

        // when using floating point samples, use a scaler instead of a divider
        // because division is much slower operation than multiplying.
        double dScaler = 1.0 / (double)resultDivider;

        Debug.Assert(length != 0);
        Debug.Assert(src.Length > 0);
        Debug.Assert(filterCoeffs.Length > 0);
        Debug.Assert(numChannels < 16);

        end = numChannels * (numSamples - length);
        outputSampleNum = numSamples - length;
        float[] dest = new float[end + numChannels - 1];

        for (j = 0; j < end; j += numChannels)
        {
            double[] sums = new double[16];
            int c, i;

            for (c = 0; c < numChannels; c++)
            {
                sums[c] = 0;
            }

            for (i = 0; i < length; i++)
            {
                float coef = filterCoeffs[i];
                for (c = 0; c < numChannels; c++)
                {
                    sums[c] += src[j + i * numChannels + c] * coef;
                }
            }

            for (c = 0; c < numChannels; c++)
            {
                sums[c] *= dScaler;
                dest[j + c] = (float)sums[c];
            }
        }
        return dest;
    }
}