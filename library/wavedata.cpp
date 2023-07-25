const char* TEXT_DATA_DOES_NOT_MATCH_DESCRIPTION = "Structure of file does not match delimited format requested by application.";
const char* TIME_OUT_OF_RANGE = "Requested time value out of range of wave.";
const char* INDEX_OUT_OF_RANGE = "Requested data index out of range of wave.";
const char* DT_MUST_BE_POSITIVE = "Delta-t must be positive, nonzero.";
const char* NO_INPUT_STREAM = "Input stream does not exist.";
const char* NO_OUTPUT_STREAM = "Output stream does not exist.";

void strcpy_pad_null(char* dest, const char* src, int destsize)
{
    bool gotnull = false;
    for(int i=0; i<destsize; i++)
    {
        if(gotnull)
            dest[i] = 0;
        else
        {
            dest[i] = src[i];
            if(src[i] == 0)
                gotnull = true;
        }
    }
}
