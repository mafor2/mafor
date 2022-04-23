# mafor
MAFOR v2 - community model for aerosol dynamics

The MAFOR (Multicomponent Aerosol FORmation) model is a Lagrangian type sectional aerosol box model which includes gas phase and aqueous phase chemistry in addition to aerosol dynamics. MAFOR simultaneously solves the time evolution of both the particle number and the mass concentrations of aerosol components in each size section in a consistent manner. The model allows for the changes in the average density of particles and represents the growth of particles in terms of both the particle number and mass. The aerosol dynamics in MAFOR are coupled to a detailed gas-phase chemistry module, which offers full flexibility for inclusion of new chemical species and reactions.

The software program is a console application with no GUI at the moment. The main intention of distributing the MAFOR model is the usage for educational and research purposes. We are looking for interested developers who want to contribute to the further development of MAFOR. All code is written in Fortran.

Model users preparing publications resulting from the usage of MAFOR are requested to cite:
    1.  Karl, M., Pirjola, L., Grönholm, T., Kurppa, M., Anand, S.,
        Zhang, X., Held, A., Sander, R., Dal Maso, M., Topping, D., 
        Jiang, S., Kangas, L., and J. Kukkonen, Description and evaluation 
        of the community aerosol dynamics model MAFOR v2.0, Geosci. Model 
        Dev. Discuss. [preprint], doi:10.5194/gmd-2021-397, 2021.
    2.  Karl, M., Kukkonen, J., Keuken, M.P., Lutzenkirchen, S.,
        Pirjola, L., Hussein, T., Modelling and measurements of urban
        aerosol processes on the neighborhood scale in Rotterdam,
        Oslo and Helsinki, Atmos. Chem. Phys., 16,
        4817-4835, doi:10.5194/acp-16-4817-2016, 2016.
