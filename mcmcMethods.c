#include "mcmcMethods.h"



struct f_params {
    double omegam; 
    double alpha; 
    double z;
};


struct int_params {
    double omegam; 
    double alpha;   
};

void progressbar(int current, int total) {
    float percent = 100 * (float)current / (float)total;
    int filled = 50 * current / total;
    //char bar[] = filled * '█' + '-' * (50 - filled);

    printf("\r%i/%i |", current, total);

    for (int i = 0; i < filled; i++)
    {   
        printf("\r");
        printf("█");
    }
    for (int i = 0; i < (50 - filled); i++)
    {
        printf("_");
    }
    
    printf("| %f%% Complete\r", percent);
    
    if (current==total) {
        printf("\n");
    }
}

void progress_simple(int current, int total) {
    printf("\r%i/%i\r", current, total);
}


void load_data(FILE *fp, double *arr) {
    
    int line = 0;

    double buffer;
     
    while (fscanf(fp, "%lf", &buffer) > 0)
    {
        arr[line] = buffer;
        line++;
    }
    
    fclose(fp);
}


double f(double x, void *params) {
    
    struct f_params *p;
    p = (struct f_params *)params;

    double omegam = (p->omegam);
    double alpha = (p->alpha);
    double z = (p->z);

    return (x * x - pow(x, alpha) * (1. - omegam) - omegam * (1. + z) * (1. + z) * (1. + z));
}


double df(double x, void *params) {

    struct f_params *p;
    p = (struct f_params *)params;    

    double omegam = (p->omegam);
    double alpha = (p->alpha);

    return (2. * x - alpha * pow(x, (alpha - 1.)) * (1. - omegam));
}

double xInit(double omegam, double z) {
    return (sqrt(omegam * pow((1. + z), 3.0) + 1. - omegam));
}


void fdf(double x, void *params, double *y, double *dy) {
    
    struct f_params *p;
    p = (struct f_params *)params;    

    double omegam = (p->omegam);
    double alpha = (p->alpha);
    double z = (p->z);
    
    double a = pow(x, alpha - 1.) * (1. - omegam);
    double zp = 1. + z;

    *y = x * x - x * a - omegam * zp * zp * zp;
    *dy = 2* x - a * alpha;
    return;
}


double integrand(double z, void *iparams) {

    struct int_params *p;
    p = (struct int_params *)iparams;    

    double omegam = (p->omegam);
    double alpha = (p->alpha);
 

    int status;
    int i;
    const gsl_root_fdfsolver_type *solver_type;
    gsl_root_fdfsolver *solver;
    gsl_function_fdf FDF;
    double x0;
    double x = xInit(omegam, z);

    struct f_params params = { omegam, alpha, z };


    FDF.f = &f;
    FDF.df = &df;
    FDF.fdf = &fdf;
    FDF.params = &params;

    solver_type = gsl_root_fdfsolver_newton;
    solver = gsl_root_fdfsolver_alloc(solver_type);
    gsl_root_fdfsolver_set(solver, &FDF, x);
    
    status = GSL_CONTINUE;
    for (i = 1; i <= MAX_ITERATIONS && status == GSL_CONTINUE; ++i) {
        /* Iterate one step of the solver */
        status = gsl_root_fdfsolver_iterate(solver);

        /* Get the new approximate solution */
        x0 = x;
        x = gsl_root_fdfsolver_root(solver);

        
        /* Check to see if the solution is within 0.00005 */
        status = gsl_root_test_delta(x, x0, 0, DREL);
        
    }
    
   

    /*
    if (status == GSL_CONTINUE) {
        printf("error: too many iterations");
    } else if (status != GSL_SUCCESS) {
        printf("error: %s\n", gsl_strerror(status));
    } */

    gsl_root_fdfsolver_free(solver);
    

    return (1 / (x));
}

double dl(double z, double omegam, double alpha, double H0) {
    
    struct int_params params = {omegam, alpha};

    gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(65536);

    double result; 
    double error;
    size_t nEvals = 40000;

    gsl_function FUNC;
    FUNC.function = &integrand;
    FUNC.params = &params;

    //gsl_set_error_handler_off();

    //gsl_integration_qng(&FUNC, 0.0, z, 0, DREL, &result, &error, &nPoints);
    gsl_integration_cquad (&FUNC, 0.0, z, 0, DREL, w, &result, &error, &nEvals);


    gsl_integration_cquad_workspace_free (w);

    //printf("check\n");

    return C0 * (1. + z) * (result / H0);
}

double mtheo(double z, double omegam, double alpha, double H0) {
    return 25. + 5 * log10(dl(z, omegam, alpha, H0));
}

double chi2(double omegam, double H0, double alpha, double *z_obs, double *mb, double *e_mb) {
    
    double sum = 0.;
    double numerator;
    for (int i = 0; i < NDAT-1; i++)
    {   
        //printf("%i\n", i);
        numerator = mtheo(z_obs[i], omegam, alpha, H0) - mb[i];
        sum += (numerator * numerator) / (e_mb[i] * e_mb[i]);
    }
    
    return sum;
}


int check(double *omegam, double *H0, double *alpha) {
    
    if (*omegam < 0) {
        return 0;   
    }    
    if (*alpha >= 2.) {
        return 0;  
    } else if (*alpha < -8)
    {
        return 0;
    }
    
    return 1;
}


double randn() {
    return (rand()+1.0)/(RAND_MAX+1.0);
}


double randnorm(double loc, double scale) {
  double y = sqrt(-2*log(randn())) * cos(2*PI*randn());
  return (loc + scale * y);
}

double min(double x, double y) {
    if (x <= y) {
        return x;
    } else {
        return y;
    }   
}

double **chain(int steps, double *z_obs, double *mb, double *e_mb) {
    
    double* values = (double*)malloc(3 * steps * sizeof(double));
    double** rows = (double**)malloc(3 * sizeof(double*));                  //creating arr[3][steps]

    for (int i=0; i<3; ++i)
    {
        rows[i] = values + i * steps;
    }
    
        double omegam = randnorm(0.3, 0.2);
        double H0 = randnorm(70., 5.);
        double alpha = randnorm(0., 1.);

    while (!check(&omegam, &H0, &alpha))
    {
        omegam = randnorm(0.3, 0.2);
        H0 = randnorm(70., 5.);
        alpha = randnorm(0., 1.);
    }
    

    double omegaOld = omegam;
    double HOld = H0;
    double alphaOld = alpha;
    
    double chiOld = chi2(omegaOld, HOld, alphaOld, z_obs, mb, e_mb);

    double chiNew, pTest;

    for (int i = 0; i < steps; i++)
    {   
        //progressbar(i, steps);
        progress_simple(i, steps);
        omegam = randnorm(omegaOld, 0.03);
        H0 = randnorm(HOld, 0.4);
        alpha = randnorm(alphaOld, 0.5);

        if (check(&omegam, &H0, &alpha)) {
            chiNew = chi2(omegam, H0, alpha, z_obs, mb, e_mb);
        } else
        {
            chiNew = 1e30;
        }
        

        pTest = exp(-0.5 * (chiNew - chiOld));
        if (min(pTest, 1.0) >= randn())
        {   
            rows[0][i] = omegam;
            rows[1][i] = H0;
            rows[2][i] = alpha;
            
            omegaOld = omegam;
            HOld = H0;
            alphaOld = alpha;

            chiOld = chiNew;
        } else
        {
            rows[0][i] = omegaOld;
            rows[1][i] = HOld;
            rows[2][i] = alphaOld;
        }
    }

    return rows;
}   

