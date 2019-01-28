//----------------------------------------------
//	CreateTime  : #CREATIONDATE#
//	Author      : Taylor Liang
//	Project     : #PROJECTNAME#
//	Company     : #SMARTDEVELOPERS#
//	Instruction : RateTransposer
//	ChangeLog   : None
//----------------------------------------------

public class RateTransposer
{
    public int Channels;

    //public ITransposer Transposer = new LinearTransposer();
    //public ITransposer Transposer = new CubicTransposer();
    public ITransposer Transposer = new ShannonTransposer();
    public AAFilter AAFilter = new AAFilter();

    private double Rate = 1;

    public void SetRate(double p_rate)
    {
        Rate = p_rate;

        Transposer.SetRate(p_rate);

        double fCutoff;
        // design a new anti-alias filter
        if (p_rate > 1.0)
        {
            fCutoff = 0.5 / p_rate;
        }
        else
        {
            fCutoff = 0.5 * p_rate;
        }
        AAFilter.SetCutoffFreq(fCutoff);
    }

    public float[] Process(float[] p_buffer, ref int outputSampleNum)
    {
        if (Rate < 1)
        {
            var transposed = Transposer.Transpose(p_buffer, Channels, p_buffer.Length / Channels, ref outputSampleNum);

            return AAFilter.Evaluate(transposed, Channels, transposed.Length / Channels, ref outputSampleNum);
        }
        else
        {
            var filtered = AAFilter.Evaluate(p_buffer, Channels, p_buffer.Length / Channels, ref outputSampleNum);

            return Transposer.Transpose(filtered, Channels, filtered.Length / Channels, ref outputSampleNum);
        }

        //return Transposer.Transpose(p_buffer, Channels, p_buffer.Length / Channels, ref outputSampleNum);
    }
}