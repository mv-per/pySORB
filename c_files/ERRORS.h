#ifndef ERRORS_H
#define ERRORS_H

/*
    Sources:
    1 - 10.1016/j.cej.2009.09.013
    2 - 10.1006/jcis.2002.8664
    3 - 10.1021/la800725s
*/


//n = number of data points
// p = number of parameters of the isotherm
double  get_error(double* n_exp, double* n_calc, size_t n, size_t p,string ERROR_FUNC){
    // printf("%f, %f \n",n_calc[0],n_exp[0])
    size_t i;
    double FOBJ = 0.;
    // double FOBJ = 0.;

    // Sum squares errors
    if (ERROR_FUNC == "SSE" || ERROR_FUNC == "ERRSQ"){
        for (i = 0; i<n;i++){
            FOBJ += pow(n_calc[i]-n_exp[i],2.0);
            // printf("%f, %f , %f\n",n_calc[i],n_exp[i], FOBJ);
        }  
    }

    // Hybrid fractional error function
    else if (ERROR_FUNC == "HYBRID"){
        for (i = 0; i<n;i++){
            FOBJ += pow(n_exp[i]-n_calc[i],2)/n_exp[i];
            // printf("%f, %f , %f\n",n_calc[i],n_exp[i], FOBJ);
        }  
        FOBJ =  100./ (double(n) -double(p)) * FOBJ;
    }

    // Average relative error
    else if (ERROR_FUNC == "ARE")
    {
        for (i = 0; i<n;i++){
            FOBJ += (n_exp[i]-n_calc[i])/n_exp[i];
        }  
        // printf("%f \t %f \n",100./ double(n), FOBJ);
        FOBJ = 100./ double(n) * FOBJ;
        
    }

    // Sum of absolute error
    else if (ERROR_FUNC == "EABS")
    {
        for (i = 0; i<n;i++)
        {
            FOBJ += fabs(n_exp[i]-n_calc[i]);
        }  
    }
    
    // Marquardt’s percent standard deviation
    else if (ERROR_FUNC == "MPSD")
    {
        for (i = 0; i<n;i++){
            FOBJ +=  pow((n_exp[i]-n_calc[i])/n_exp[i],2.0);
        }  
        FOBJ = 100. * sqrt(1./ (double(n) -double(p)) * FOBJ);
    }

    // Standard deviation of relative errors
    
    else if (ERROR_FUNC == "SDRE") 
    {
        for (i = 0; i<n;i++){
            FOBJ +=  pow((n_exp[i]-n_calc[i])/n_exp[i] - get_error(n_exp, n_calc, n, p,"ARE"),2.0);
        }  
        FOBJ = sqrt(FOBJ/(double(n)-1.));
    }

    // Spearman's correlation coefficient
    else if (ERROR_FUNC == "R_S") 
    {
        FOBJ = 1 - 6.0 * get_error(n_exp, n_calc, n, p,"SSE") / double (n) / pow(double(n)-1.0, 2.0);
    }
    
    
    // Nonlinear chi-square test
    else if (ERROR_FUNC == "CHI_2")
    {
        for (i = 0; i<n;i++){
            FOBJ += pow(n_exp[i]-n_calc[i],2.0)/n_exp[i];
        }  
    }
  
    else{
        printf("Error Function not found, returning Average relative error");
        for (i = 0; i<n;i++){
            FOBJ += (n_exp[i]-n_calc[i])/n_exp[i];
        }  
        FOBJ =  100./double (n) * FOBJ;
    }
    // printf("%f \n", FOBJ);
    // if (isinf(FOBJ) || isnan(FOBJ)){
    //     return 1000.;
    // }
    // else {return FOBJ;}

    return FOBJ;
}

#endif