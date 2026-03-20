This is the yuckiest bit of hacking I have ever done, but it works. I first tried to patch individual versions of SCILAB, hoping that the SCILAB developers would add them to later versions. Then came Version 6. The thread handling for the V6 parser was crap and it completely stuffed the Event Handler routines I used for mouse interactions. So I just gave up and fixed LaueG to only work with the final Version 5 of SCILAB, v5.5.2.

There is no need for you to recompile these modules, but I include the details here mainly so you can feel my pain.
