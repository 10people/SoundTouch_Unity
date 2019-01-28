//----------------------------------------------
//	CreateTime  : #CREATIONDATE#
//	Author      : Taylor Liang
//	Project     : #PROJECTNAME#
//	Company     : #SMARTDEVELOPERS#
//	Instruction : TDStretch
//	ChangeLog   : None
//----------------------------------------------

using System;
using UnityEngine;

public class TDStretch
{
    public double Tempo = 1;
    public int Channels;
    public int SampleRate;

    private int seekWindowLength;
    private int seekLength;
    private int seekWindowMs;
    private int sequenceMs;
    private int overlapLength;

    private double nominalSkip;
    private int sampleReq;
    private double skipFract;

    public float[] Process(float[] inputBuffer)
    {
        skipFract = 0;

        int ovlSkip;
        int offset = 0;
        int temp;
        float[] outputBuffer = new float[0];
        int inputBufferIndex = 0;
        float[] compareBuffer = new float[Channels * overlapLength];

        bool isFirstTime = true;

        while (inputBuffer.Length - inputBufferIndex >= sampleReq)
        {
            if (isFirstTime)
            {
                //first time
                //int skip = 0;
                int skip = (int)(Tempo * overlapLength + 0.5);
                skipFract -= skip;
                Debug.Assert(nominalSkip >= -skipFract);

                isFirstTime = false;
            }
            else
            {

                // apart from the very beginning of the track, 
                // scan for the best overlapping position & do overlap-add
                offset = SeekBestOverlapPosition(inputBuffer, compareBuffer, inputBufferIndex);

                // Mix the samples in the 'inputBuffer' at position of 'offset' with the 
                // samples in 'midBuffer' using sliding overlapping
                // ... first partially overlap with the end of the previous sequence
                // (that's in 'midBuffer')
                var result = Overlap(inputBuffer, compareBuffer, offset);
                AddBuffer(ref outputBuffer, result, 0, result.Length);
                offset += overlapLength;
            }

            if (inputBuffer.Length - inputBufferIndex >= (offset + seekWindowLength - overlapLength))
            {
                // length of sequence
                temp = (seekWindowLength - 2 * overlapLength);

                AddBuffer(ref outputBuffer, inputBuffer, inputBufferIndex + Channels * offset, temp);

                // Copies the end of the current sequence from 'inputBuffer' to 
                // 'midBuffer' for being mixed with the beginning of the next 
                // processing sequence and so on
                Debug.Assert((offset + temp + overlapLength) <= inputBuffer.Length - inputBufferIndex);
                Array.Copy(inputBuffer, inputBufferIndex + Channels * (offset + temp), compareBuffer, 0, Channels * overlapLength);

                // Remove the processed samples from the input buffer. Update
                // the difference between integer & nominal skip step to 'skipFract'
                // in order to prevent the error from accumulating over time.
                skipFract += nominalSkip;   // real skip size
                ovlSkip = (int)skipFract;   // rounded to integer skip
                skipFract -= ovlSkip;       // maintain the fraction part, i.e. real vs. integer skip
                inputBufferIndex += ovlSkip;
            }
        }

        return outputBuffer;
    }

    float[] Overlap(float[] inputBuffer, float[] compareBuffer, int inputOffset)
    {
        if (Channels == 1)
        {
            // mono sound.
            return OverlapMono(inputBuffer, compareBuffer, inputOffset);
        }
        else if (Channels == 2)
        {
            // stereo sound
            return OverlapStereo(inputBuffer, compareBuffer, 2 * inputOffset);
        }
        else
        {
            Debug.Assert(Channels > 0);
            return OverlapMulti(inputBuffer, compareBuffer, Channels * inputOffset);
        }
    }

    float[] OverlapMono(float[] inputBuffer, float[] compareBuffer, int inputOffset)
    {
        int i;
        float m1, m2;
        float[] result = new float[overlapLength];

        m1 = 0;
        m2 = overlapLength;

        for (i = 0; i < overlapLength; i++)
        {
            result[i] = (inputBuffer[inputOffset + i] * m1 + compareBuffer[i] * m2) / overlapLength;
            m1 += 1;
            m2 -= 1;
        }

        return result;
    }

    // Overlaps samples in 'midBuffer' with the samples in 'pInput'
    float[] OverlapStereo(float[] inputBuffer, float[] compareBuffer, int inputOffset)
    {
        int i;
        float fScale;
        float f1;
        float f2;
        float[] result = new float[2 * overlapLength];

        fScale = 1.0f / (float)overlapLength;

        f1 = 0;
        f2 = 1.0f;

        for (i = 0; i < 2 * (int)overlapLength; i += 2)
        {
            result[i + 0] = inputBuffer[inputOffset + i + 0] * f1 + compareBuffer[i + 0] * f2;
            result[i + 1] = inputBuffer[inputOffset + i + 1] * f1 + compareBuffer[i + 1] * f2;

            f1 += fScale;
            f2 -= fScale;
        }

        return result;
    }


    // Overlaps samples in 'midBuffer' with the samples in 'input'. 
    float[] OverlapMulti(float[] inputBuffer, float[] compareBuffer, int inputOffset)
    {
        int i;
        float fScale;
        float f1;
        float f2;
        float[] result = new float[2 * overlapLength];

        fScale = 1.0f / (float)overlapLength;

        f1 = 0;
        f2 = 1.0f;

        i = 0;
        for (int i2 = 0; i2 < overlapLength; i2++)
        {
            // note: Could optimize this slightly by taking into account that always channels > 2
            for (int c = 0; c < Channels; c++)
            {
                result[i] = inputBuffer[inputOffset + i] * f1 + compareBuffer[i] * f2;
                i++;
            }
            f1 += fScale;
            f2 -= fScale;
        }

        return result;
    }

    void AddBuffer(ref float[] output, float[] input, int inputIndex, int copyNum)
    {
        var outputIndex = output.Length;
        Array.Resize(ref output, output.Length + copyNum);
        Array.Copy(input, inputIndex, output, outputIndex, copyNum);
    }

    private const int SCANSTEP = 16;
    private const int SCANWIND = 8;
    private const float FLT_MAX = 3.402823466e+38F;     // max value

    public int SeekBestOverlapPosition(float[] inputBuffer, float[] compareBuffer, int seekStart)
    {
        //return SeekBestOverlapPositionQuick(inputBuffer, compareBuffer, seekStart);
        return SeekBestOverlapPositionFull(inputBuffer, compareBuffer, seekStart);
    }

    private int SeekBestOverlapPositionFull(float[] inputBuffer, float[] compareBuffer, int seekStart)
    {
        int bestOffs;
        double bestCorr;
        int i;
        double norm = 0;

        bestCorr = -FLT_MAX;
        bestOffs = 0;

        // Scans for the best correlation value by testing each possible position
        // over the permitted range.
        bestCorr = calcCrossCorr(inputBuffer, seekStart, compareBuffer, ref norm);
        bestCorr = (bestCorr + 0.1) * 0.75;

        for (i = 1; i < seekLength; i++)
        {
            double corr;

            // Calculates correlation value for the mixing position corresponding to 'i'
            // In non-parallel version call "calcCrossCorrAccumulate" that is otherwise same
            // as "calcCrossCorr", but saves time by reusing & updating previously stored 
            // "norm" value
            corr = calcCrossCorrAccumulate(inputBuffer, seekStart + Channels * i, compareBuffer, ref norm);

            // heuristic rule to slightly favour values close to mid of the range
            double tmp = (double)(2 * i - seekLength) / (double)seekLength;
            corr = ((corr + 0.1) * (1.0 - 0.25 * tmp * tmp));

            // Checks for the highest correlation value
            if (corr > bestCorr)
            {
                // For optimal performance, enter critical section only in case that best value found.
                // in such case repeat 'if' condition as it's possible that parallel execution may have
                // updated the bestCorr value in the mean time
                if (corr > bestCorr)
                {
                    bestCorr = corr;
                    bestOffs = i;
                }
            }
        }

        // clear cross correlation routine state if necessary (is so e.g. in MMX routines).
        //default is empty

        return bestOffs;
    }

    private int SeekBestOverlapPositionQuick(float[] inputBuffer, float[] compareBuffer, int seekStart)
    {
        int bestOffs;
        int i;
        int bestOffs2;
        float bestCorr, corr;
        float bestCorr2;
        double norm = 0;

        // note: 'float' types used in this function in case that the platform would need to use software-fp

        bestCorr = bestCorr2 = -FLT_MAX;
        bestOffs = bestOffs2 = SCANWIND;

        // Scans for the best correlation value by testing each possible position
        // over the permitted range. Look for two best matches on the first pass to
        // increase possibility of ideal match.
        //
        // Begin from "SCANSTEP" instead of SCANWIND to make the calculation
        // catch the 'middlepoint' of seekLength vector as that's the a-priori 
        // expected best match position
        //
        // Roughly:
        // - 15% of cases find best result directly on the first round,
        // - 75% cases find better match on 2nd round around the best match from 1st round
        // - 10% cases find better match on 2nd round around the 2nd-best-match from 1st round
        for (i = SCANSTEP; i < seekLength - SCANWIND - 1; i += SCANSTEP)
        {
            // Calculates correlation value for the mixing position corresponding
            // to 'i'
            corr = (float)calcCrossCorr(inputBuffer, seekStart + Channels * i, compareBuffer, ref norm);
            // heuristic rule to slightly favour values close to mid of the seek range
            float tmp = (float)(2 * i - seekLength - 1) / (float)seekLength;
            corr = ((corr + 0.1f) * (1.0f - 0.25f * tmp * tmp));

            // Checks for the highest correlation value
            if (corr > bestCorr)
            {
                // found new best match. keep the previous best as 2nd best match
                bestCorr2 = bestCorr;
                bestOffs2 = bestOffs;
                bestCorr = corr;
                bestOffs = i;
            }
            else if (corr > bestCorr2)
            {
                // not new best, but still new 2nd best match
                bestCorr2 = corr;
                bestOffs2 = i;
            }
        }

        // Scans surroundings of the found best match with small stepping
        int end = Math.Min(bestOffs + SCANWIND + 1, seekLength);
        for (i = bestOffs - SCANWIND; i < end; i++)
        {
            if (i == bestOffs) continue;    // this offset already calculated, thus skip

            // Calculates correlation value for the mixing position corresponding
            // to 'i'
            corr = (float)calcCrossCorr(inputBuffer, seekStart + Channels * i, compareBuffer, ref norm);
            // heuristic rule to slightly favour values close to mid of the range
            float tmp = (float)(2 * i - seekLength - 1) / (float)seekLength;
            corr = ((corr + 0.1f) * (1.0f - 0.25f * tmp * tmp));

            // Checks for the highest correlation value
            if (corr > bestCorr)
            {
                bestCorr = corr;
                bestOffs = i;
            }
        }

        // Scans surroundings of the 2nd best match with small stepping
        end = Math.Min(bestOffs2 + SCANWIND + 1, seekLength);
        for (i = bestOffs2 - SCANWIND; i < end; i++)
        {
            if (i == bestOffs2) continue;    // this offset already calculated, thus skip

            // Calculates correlation value for the mixing position corresponding
            // to 'i'
            corr = (float)calcCrossCorr(inputBuffer, seekStart + Channels * i, compareBuffer, ref norm);
            // heuristic rule to slightly favour values close to mid of the range
            float tmp = (float)(2 * i - seekLength - 1) / (float)seekLength;
            corr = ((corr + 0.1f) * (1.0f - 0.25f * tmp * tmp));

            // Checks for the highest correlation value
            if (corr > bestCorr)
            {
                bestCorr = corr;
                bestOffs = i;
            }
        }

        return bestOffs;
    }

    private double calcCrossCorrAccumulate(float[] input, int inputOffsetIndex, float[] compare, ref double norm)
    {
        double corr;
        int i;

        corr = 0;

        // cancel first normalizer tap from previous round
        for (i = 1; i <= Channels; i++)
        {
            norm -= input[-i + inputOffsetIndex] * input[-i + inputOffsetIndex];
        }

        // Same routine for stereo and mono. For Stereo, unroll by factor of 2.
        // For mono it's same routine yet unrollsd by factor of 4.
        for (i = 0; i < Channels * overlapLength; i += 4)
        {
            corr += input[i + inputOffsetIndex] * compare[i] +
                    input[i + inputOffsetIndex + 1] * compare[i + 1] +
                    input[i + inputOffsetIndex + 2] * compare[i + 2] +
                    input[i + inputOffsetIndex + 3] * compare[i + 3];
        }

        // update normalizer with last samples of this round
        for (int j = 0; j < Channels; j++)
        {
            i--;
            norm += input[i + inputOffsetIndex] * input[i + inputOffsetIndex];
        }

        return corr / Math.Sqrt((norm < 1e-9 ? 1.0 : norm));
    }

    private double calcCrossCorr(float[] input, int inputOffsetIndex, float[] compare, ref double anorm)
    {
        double corr;
        double norm;
        int i;

        corr = norm = 0;
        // Same routine for stereo and mono. For Stereo, unroll by factor of 2.
        // For mono it's same routine yet unrollsd by factor of 4.
        for (i = 0; i < Channels * overlapLength; i += 4)
        {
            corr += input[i + inputOffsetIndex] * compare[i] +
                    input[i + inputOffsetIndex + 1] * compare[i + 1];

            norm += input[i + inputOffsetIndex] * input[i + inputOffsetIndex] +
                    input[i + inputOffsetIndex + 1] * input[i + inputOffsetIndex + 1];

            // unroll the loop for better CPU efficiency:
            corr += input[i + inputOffsetIndex + 2] * compare[i + 2] +
                    input[i + inputOffsetIndex + 3] * compare[i + 3];

            norm += input[i + inputOffsetIndex + 2] * input[i + inputOffsetIndex + 2] +
                    input[i + inputOffsetIndex + 3] * input[i + inputOffsetIndex + 3];
        }

        anorm = norm;
        return corr / Math.Sqrt((norm < 1e-9 ? 1.0 : norm));
    }

    // Adjust tempo param according to tempo, so that variating processing sequence length is used
    // at various tempo settings, between the given low...top limits
    private const float AUTOSEQ_TEMPO_LOW = 0.5f;    // auto setting low tempo range (-50%)
    private const float AUTOSEQ_TEMPO_TOP = 2.0f;   // auto setting top tempo range (+100%)

    // sequence-ms setting values at above low & top tempo
    private const float AUTOSEQ_AT_MIN = 90.0f;
    private const float AUTOSEQ_AT_MAX = 40.0f;

    // seek-window-ms setting values at above low & top tempoq
    private const float AUTOSEEK_AT_MIN = 20.0f;
    private const float AUTOSEEK_AT_MAX = 15.0f;
    public void CalcSeqParameters()
    {
        double seq, seek;

        //bAutoSeqSetting)
        seq = (AUTOSEQ_AT_MIN - (((AUTOSEQ_AT_MAX - AUTOSEQ_AT_MIN) / (AUTOSEQ_TEMPO_TOP - AUTOSEQ_TEMPO_LOW))) * (AUTOSEQ_TEMPO_LOW)) + ((AUTOSEQ_AT_MAX - AUTOSEQ_AT_MIN) / (AUTOSEQ_TEMPO_TOP - AUTOSEQ_TEMPO_LOW)) * Tempo;
        seq = Math_Clamp(seq, AUTOSEQ_AT_MAX, AUTOSEQ_AT_MIN);
        sequenceMs = (int)(seq + 0.5);

        //bAutoSeekSetting)
        seek = (AUTOSEEK_AT_MIN - (((AUTOSEEK_AT_MAX - AUTOSEEK_AT_MIN) / (AUTOSEQ_TEMPO_TOP - AUTOSEQ_TEMPO_LOW))) * (AUTOSEQ_TEMPO_LOW)) + ((AUTOSEEK_AT_MAX - AUTOSEEK_AT_MIN) / (AUTOSEQ_TEMPO_TOP - AUTOSEQ_TEMPO_LOW)) * Tempo;
        seek = Math_Clamp(seek, AUTOSEEK_AT_MAX, AUTOSEEK_AT_MIN);
        seekWindowMs = (int)(seek + 0.5);

        // Update seek window lengths
        seekWindowLength = (SampleRate * sequenceMs) / 1000;
        if (seekWindowLength < 2 * overlapLength)
        {
            seekWindowLength = 2 * overlapLength;
        }
        seekLength = (SampleRate * seekWindowMs) / 1000;
    }

    public void CalcTempoRelated(double p_tempo)
    {
        int intskip;

        Tempo = p_tempo;

        // Calculate new sequence duration
        CalcSeqParameters();

        // Calculate ideal skip length (according to tempo value) 
        nominalSkip = Tempo * (seekWindowLength - overlapLength);
        intskip = (int)(nominalSkip + 0.5);

        // Calculate how many samples are needed in the 'inputBuffer' to 
        // process another batch of samples
        //sampleReq = max(intskip + overlapLength, seekWindowLength) + seekLength / 2;
        sampleReq = Math.Max(intskip + overlapLength, seekWindowLength) + seekLength;
    }

    public void SetOverlapLength(int aOverlapMS)
    {
        CalculateOverlapLength(aOverlapMS);
    }

    void CalculateOverlapLength(int overlapInMsec)
    {
        int newOvl;

        newOvl = (SampleRate * overlapInMsec) / 1000;
        if (newOvl < 16) newOvl = 16;

        // must be divisible by 8
        newOvl -= newOvl % 8;

        overlapLength = newOvl;
    }

    private double Math_Clamp(double x, double a, double b)
    {
        return Math.Min(Math.Max(x, a), b);
    }
}