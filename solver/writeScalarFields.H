bool writeTime = false;

if (word(runTime.controlDict().lookup("writeControl")) == "timeStep")
{
    writeTime =
   !(
        runTime.timeIndex()
       %readLabel(runTime.controlDict().lookup("writeInterval"))
    );
}

if (word(runTime.controlDict().lookup("writeControl")) == "runTime")
{
    label writeIndex = label
    (
        (
            (
                runTime.timeOutputValue()
              - runTime.startTime().value()
            )
          + 0.5*runTime.deltaT().value()
        )
       /readLabel(runTime.controlDict().lookup("writeInterval"))
    );

    if (writeIndex > writeTimeIndex)
    {
        writeTime = true;
        writeTimeIndex = writeIndex;
    }
}

if (writeTime)
{
    os  << runTime.timeOutputValue() << "\t"
        << y[0] << "\t"
        << y[1] << endl;
}
