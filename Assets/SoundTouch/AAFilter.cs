//----------------------------------------------
//	CreateTime  : 14/01/2019 17:37:50
//	Author      : Taylor Liang
//	Project     : My_Demo
//	Company     : BuGu
//	Instruction : AAFilter
//	ChangeLog   : None
//----------------------------------------------

using System;
using UnityEngine;

public class AAFilter
{
    public AAFilter()
    {
        CalculateCoeffs();
    }

    public FIRFilter FirFilter = new FIRFilter();

    /// Low-pass filter cut-off frequency, negative = invalid
    double cutoffFreq;

    /// num of filter taps
    int length = 64;

    public float[] Evaluate(float[] src, int channels, int sampleNum, ref int outputSampleNum)
    {
        return FirFilter.Evaluate(src, sampleNum, channels, ref outputSampleNum);
    }

    // Calculates coefficients for a low-pass FIR filter using Hamming window
    void CalculateCoeffs()
    {
        int i;
        double cntTemp, temp, tempCoeff, h, w;
        double wc;
        double scaleCoeff, sum;
        double[] work = new double[length];
        float[] coeffs = new float[length];

        Debug.Assert(length >= 2);
        Debug.Assert(length % 4 == 0);
        Debug.Assert(cutoffFreq >= 0);
        Debug.Assert(cutoffFreq <= 0.5);

        wc = 2.0 * Math.PI * cutoffFreq;
        tempCoeff = Math.PI * 2 / (double)length;

        sum = 0;
        for (i = 0; i < length; i++)
        {
            cntTemp = (double)i - (double)(length / 2);

            temp = cntTemp * wc;
            if (temp != 0)
            {
                h = Math.Sin(temp) / temp;                     // sinc function
            }
            else
            {
                h = 1.0;
            }
            w = 0.54 + 0.46 * Math.Cos(tempCoeff * cntTemp);       // hamming window

            temp = w * h;
            work[i] = temp;

            // calc net sum of coefficients 
            sum += temp;
        }

        // ensure the sum of coefficients is larger than zero
        Debug.Assert(sum > 0);

        // ensure we've really designed a lowpass filter...
        Debug.Assert(work[length / 2] > 0);
        Debug.Assert(work[length / 2 + 1] > -1e-6);
        Debug.Assert(work[length / 2 - 1] > -1e-6);

        // Calculate a scaling coefficient in such a way that the result can be
        // divided by 16384
        scaleCoeff = 16384.0f / sum;

        for (i = 0; i < length; i++)
        {
            temp = work[i] * scaleCoeff;
            // scale & round to nearest integer
            temp += (temp >= 0) ? 0.5 : -0.5;
            // ensure no overfloods
            Debug.Assert(temp >= -32768 && temp <= 32767);
            coeffs[i] = (float)temp;
        }

        // Set coefficients. Use divide factor 14 => divide result by 2^14 = 16384
        FirFilter.SetCoefficients(coeffs, length, 14);
    }

    public void SetCutoffFreq(double newCutoffFreq)
    {
        cutoffFreq = newCutoffFreq;
        CalculateCoeffs();
    }
}