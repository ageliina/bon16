#ifndef BON16_H
#define BON16_H

#define N_PARAMETERS 16

struct args_f_star { double log_M_star_star; double k_log_M_star_star; double alpha; double z1; };
struct args_f_lambda_SAR { double log_lambda_SAR_star0; double k_lambda; double log_M_star0; double gamma10; double k_gamma; double k_gamma1; double gamma2; double z0; double z1; };
struct args_f_z { double p1; double p2; double p3; double z0; double z1; };
struct args_Psi { double log_Psi_star; struct args_f_lambda_SAR args_l; struct args_f_star args_s; struct args_f_z args_z; };
struct args_Phi_star { double log_M_star; double z; struct args_Psi args; };
struct args_Phi_SAR { double log_lambda_SAR; double z; struct args_Psi args; };
struct args_Phi_X { double log_L_X_min; double log_L_X_max; double z; struct args_Psi args; };
struct args_N_above_S { double log_F_X_lo; double log_F_X_hi; double Gamma; struct args_Psi args; };

double get_log_f_star(double log_M_star, double z, struct args_f_star args);
double get_log_f_lambda_SAR(double log_lambda_SAR, double log_M_star, double z, struct args_f_lambda_SAR args);
double get_log_f_z(double z, struct args_f_z args);
double get_log_Psi(double log_M_star, double log_lambda_SAR, double z, struct args_Psi args);
double get_Phi_star(double x[], size_t dim, void *p);
double get_Phi_SAR(double x[], size_t dim, void *p);
double get_Phi_X(double x[], size_t dim, void *p);
double get_N_above_S(double x[], size_t dim, void *p);
double get_N(double x[], size_t dim, void *p);

#endif
