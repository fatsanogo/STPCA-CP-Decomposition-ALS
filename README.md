# STPCA-CPD-ALS initialization
This repo contains the code for the full procedure initializing the CP factors using the STPCA components as describe in the paper "A Spatio-temporal CP decomposition analysis of New England region in the US" by Fatoumata Sanogo. The data can be downloaded as describe in the paper from the NCAR Climate Data Gateway website (https://doi.org/10.5065/D6SJ1JCH).

The function stpca_to_cp_init in cp_init.R is the master function that performs our full procedure.
The code tO find the spatio temporal compenent is publicly available from the author website for the paper: "Spatio-temporal principal component analysis" by Mirosław Krzyśko,Peter NijkampORCID Icon,Waldemar Ratajczak,Waldemar Wołyński &Beata Wenerska.
The fucntion STPCA.R and STPCA_RUN.R are the codes to find the spatio-temporal component. These are from the author website.
