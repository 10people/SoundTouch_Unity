//----------------------------------------------
//	CreateTime  : #CREATIONDATE#
//	Author      : Taylor Liang
//	Project     : #PROJECTNAME#
//	Company     : #SMARTDEVELOPERS#
//	Instruction : ITransposer
//	ChangeLog   : None
//----------------------------------------------

public interface ITransposer
{
    void SetRate(double p_rate);

    float[] Transpose(float[] inputBuffer, int channels, int sampleNum, ref int outputSampleNum);
}