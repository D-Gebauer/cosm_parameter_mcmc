#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_hyperg.h>
#include <time.h>

#define MAX_ITERATIONS 5000
#define DREL 0.00005
#define DABS 0.000001
#define C0 299792.458
#define PI 3.14159265
#define NDAT 580
#define STEPS 10000


struct f_params;

struct int_params;

void progressbar(int, int);

void progress_simple(int, int);

void load_data(FILE *fp, double *arr);

double f(double x, void *params);

void fdf(double x, void *params, double *y, double *dy);

double integrand(double z, void *iparams);

double dl(double z, double omegam, double alpha, double H0);

double mtheo(double z, double omegam, double alpha, double H0);

double chi2(double omegam, double H0, double alpha, double *z_obs, double *mb, double *e_mb);

int check(double *omegam, double *H0, double *alpha);

double randn();

double randnorm(double loc, double scale);

double min(double x, double y);

double **chain(int steps, double *z_obs, double *mb, double *e_mb);
