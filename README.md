# SLEIPNIR

We propose to expand the available measurement capacity for neutron users by constructing a series of small, compact instruments at the ESS in Lund, without restricting ESS's future instrument portfolio. This would utilize beam ports that, due to spatial constraints, are unsuitable for larger state-of-the-art instruments. Up to four instruments could be installed in a limited area near a single radiation channel, each at a significantly lower cost than current ESS instrument projects. These modular instruments would not be on par with the state-of-the-art instruments being built at the ESS, but would be on par with other instruments around Europe, simply due to the high brilliance of the ESS moderator. Distinctly, not all studies require the advanced capabilities of the ESS's flagship instruments. Ambitious or high-risk experiments often need multiple preparatory trials, where the sample environment changes over a longer time frame than the measurement duration. For industrial applications, rapid access to measurement time is often more valuable than delayed access to top-tier instruments.

In this study, Monte Carlo simulation tools are used to design and optimize a version of SLEIPNIR (SLender workhorsE Instruments Pledged as Neutron Infrastructure for Research), as a proof of concept for utilizing the "dead-space" at the ESS. McStas will be used to iteratively simulate different configurations of one part of the instrument, and the design will then be bench marked against other instruments of the same type.

## Contents

### Optimization

- [SLEIPNIR PD McStas instrument](Optimization/SLEIPNER.instr)
- [Python Optimization Script](Optimization/optimize_4.py)

### Results

- Three possible versions of one leg of SLEIPNER have been optimized for different purposes. The resolutions functions are presented compared to the resolution function of D1B at the ILL in the following three figures:
    - [High Resolution](Results/D1B_4_Low.png) with a flux of 4.48e5 neutrons pr. s. pr. cm^2
    - [Medium Resolution](Results/D1B_4_Medium.png) with a flux of 2.25236e+06 neutrons pr. s. pr. cm^2
    - [Low Resolution](Results/D1B_4_High.png) with a flux of 2.76781e+07 neutrons pr. s. pr. cm^2
- Each instrument would be able to run at the designed 2.6 Å, but would also be able to acces half wavelengths with the same resolution, but approximately 1/5 of the flux.
- These modular instruments would be able to be built in parralel to maximize the use of the available space and beam time.

If you have any questions, please contact [Nicolai Amin](mailto:s194113@dtu.dk)
