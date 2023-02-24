## Known issues and work-around
### Freshwater injection stuck near injection place
  This can happen with the combined use of (1) point source/sink (`if_source`=1); (2) LSC2 with 1-2 layers near 
injection places; (3) baroclinic model. The main symptom is that the freshwater seemingly gets stuck near the 
 injection points and does not flow out as expected.
 
 The reason is that insufficient number of vertical layers cannot properly support exchange flow (that 
requires stratification), and as a result, the fresh/salt water interface oscillates instead of tilting as expects.
Some work-arounds are:

1. Use more vertical layers near injection (one way to do this is to deepen the local depths to allow more layers)
2. Nudge (inu_SAL=1 or 2) strongly in a region near injection
3. Change to open boundary condition approach
4. If the salt intrusion should never reach the injection place, change salinity initial condition by creating
   a freswater zone from injection to somewhere downstream near the max intrusion.
