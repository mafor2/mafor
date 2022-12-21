# mafor
MAFOR v2 - community model for aerosol dynamics

The MAFOR (Multicomponent Aerosol FORmation) model is a Lagrangian type sectional aerosol box model which includes gas phase and aqueous phase chemistry in addition to aerosol dynamics. The model simultaneously solves the time evolution of both the particle number and the mass concentrations of aerosol components in each size section. The model allows for the changes in the average density of particles and represents the growth of particles in terms of both the particle number and mass. The aerosol dynamics in MAFOR are coupled to a detailed gas-phase chemistry module, which offers full flexibility for inclusion of new chemical species and reactions. The MAFOR model version 2 is well documented and versatile to use, providing a range of alternative parameterizations for various aerosol processes. The model includes an efficient numerical integration of particle number and mass concentrations, an operator splitting of processes, and the use of a fixed sectional method. The model could be used as a module in various atmospheric and climatic models.

The software program is a console application with no GUI at the moment. The main intention of distributing the MAFOR model is the usage for educational and research purposes. We are looking for interested developers who want to contribute to the further development of MAFOR. All code is written in Fortran.

Model users preparing publications resulting from the usage of MAFOR are requested to cite:

    1.  Karl, M., Pirjola, L., Grönholm, T., Kurppa, M., Anand, S.,
        Zhang, X., Held, A., Sander, R., Dal Maso, M., Topping, D., 
        Jiang, S., Kangas, L., and J. Kukkonen, Description and evaluation 
        of the community aerosol dynamics model MAFOR v2.0, Geosci. Model Dev.,
        15, 3969-4026, https://doi.org/10.5194/gmd-15-3969-2022, 2022.

