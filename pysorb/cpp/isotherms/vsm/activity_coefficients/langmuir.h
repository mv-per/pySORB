// /*
//     Langmuir equations to estimate adsorbed quantities
//     based on: {b}, {n_max}, and {P}
// */

// typedef struct
// {
//     vector<vector<double>> param;
//     vector<double> y;
//     double P;
// } mix_vsm_params;

// /*
//     Returns single n_calc using Langmuir
// */
// double CalculateLangmuir(double *P, double &n_i_max, double &b)
// {
//     return n_i_max * b * *P / (1.0 + b * *P);
// }

// /*
//     Returns partial n_calc and total n_calc using Extended Langmuir
// */
// double *CalculateExtendedLangmuir(vector<vector<double>> &param, vector<double> &y, double *P, int *ncomp)
// {

// }