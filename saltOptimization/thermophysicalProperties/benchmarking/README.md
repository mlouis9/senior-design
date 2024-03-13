# Thermophysical Property Calculation Benchmarking
This directory shows the benchmarking of thermophysical property calcualtions using the custom python module: thermophysicalProperties.

## Test Cases

### Chloride RK Comparisons
1. Calculating the density of the KCl-NaCl binary system using the RK espansion. The  results are compared with the relevant plots from [[1]](#agca-empirical).
2. Calculating the density of the ternary KCl-NaCl-UCl3 system along the slice NaCl : UCl3 = 0.25 : 0.75, again the results are compared to corresponding plots from [[1]](#agca-empirical)

### Fluoride RK Comparisons
1. Calculating the density of the ternary NaF-LiF-ZrF4 system along the slice <br> NaF : LiF = 0.4 : 0.6. These results are compared to plots from [[2]](#birri-application)

## References

1. <a name="agca-empirical"></a> [Can Agca, "Empirical estimation of densities in NaCl-KCl-UCl3 and NaCl-KCl-YCl3 molten salts using Redlich-Kister expansion", Chemical Engineering Science](https://www.sciencedirect.com/science/article/pii/S0009250921006515?via%3Dihub)
2. <a name="birri-application"></a> [Anthony Birri, "Application of the Redlich-Kister expansion for estimating the density of molten fluoride psuedo-ternary salt systems of nuclear industry interest", Chemical Engineering Science](https://www.sciencedirect.com/science/article/pii/S0009250922005383#f0010)