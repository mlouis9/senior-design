# Thermophysical Property Calculation Benchmarking
This directory shows the benchmarking of thermophysical property calcualtions using the custom python module: [thermophysicalProperties](../../../modules/thermophysicalProperties.py).

## Test Cases

### Chloride RK Comparisons
1. Calculating the density of the KCl-NaCl binary system using the RK espansion. The  results are compared with the relevant plots from [[1]](#agca-empirical).
2. Calculating the density of the ternary KCl-NaCl-UCl3 system along the slice NaCl : UCl3 = 0.25 : 0.75, again the results are compared to corresponding plots from [[1]](#agca-empirical)

### Fluoride RK Comparisons
1. Calculating the density of the ternary NaF-LiF-ZrF4 system along the slice NaF : LiF = 0.4 : 0.6. These results are compared to plots from [[2]](#birri-application)

### Heat Capacity
1. The heat capacity of the following pseudo binary systems was calculated and compared to the data given in [[3]](#beilmann-excess):
   1. LiF-KF
   2. LiF-NaF

### Thermal Conductivity
1. The thermal conductivity of the KCl-MgCl$_2$ system was calculated, and compared with the corresponding plot in [[4]](#xu-experimental)

### Viscosity
1. The dynamic viscosity of LiF-NaF-KF (FLiNaK), LiF-BeF2 (66-34 mol%), and LiF-BeF2 (73-27 mol%) were calculated and compared to the reference [[5]](#tkacheva-viscosity)

## References

1. <a name="agca-empirical"></a> [Can Agca, "Empirical estimation of densities in NaCl-KCl-UCl3 and NaCl-KCl-YCl3 molten salts using Redlich-Kister expansion", Chemical Engineering Science](https://www.sciencedirect.com/science/article/pii/S0009250921006515?via%3Dihub)
2. <a name="birri-application"></a> [Anthony Birri, "Application of the Redlich-Kister expansion for estimating the density of molten fluoride psuedo-ternary salt systems of nuclear industry interest", Chemical Engineering Science](https://www.sciencedirect.com/science/article/pii/S0009250922005383#f0010)
3. <a name="beilmann-excess"></a> [M. Beilmann, "Excess Heat Capacity in Liquid Binary Alkali-Fluoride Mixtures", Inorg. Chem](https://pubs.acs.org/doi/full/10.1021/ic302168g)
4. <a name="xu-experimental"></a> [Xiankun Xu, "Experimental Test of Properties of KCl–MgCl2 Eutectic Molten Salt for Heat Transfer and Thermal Storage Fluid in Concentrated Solar Power Systems", Journal of Solar Energy Engineering](https://asmedigitalcollection.asme.org/solarenergyengineering/article/140/5/051011/368304/Experimental-Test-of-Properties-of-KCl-MgCl2)
5. <a name="tkacheva-viscosity"></a> [Olga Tkacheva, "Viscosity of fluoride melts promising for molten salt nuclear reactors", Electrochemical Materials and Technologies](https://www.researchgate.net/publication/375712728_Viscosity_of_fluoride_melts_promising_for_molten_salt_nuclear_reactors)