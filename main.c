#include "mcmcMethods.h"


int main(void) {
    
    srand(time(0));

    double z_obs[NDAT];
    double mb[NDAT];
    double e_mb[NDAT];


    load_data(fopen("./data/zUnion.txt", "r"), z_obs);
    load_data(fopen("./data/mbUnion.txt", "r"), mb);
    load_data(fopen("./data/emBUnion.txt", "r"), e_mb);


    clock_t begin = clock();
    
    double **arr = chain(STEPS, z_obs, mb, e_mb);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("time spent: %lfs\n", time_spent);


    FILE *ofp = fopen("./out/chain.txt", "w");

       
    for (int x = 0; x < STEPS; x++)
    {
        fprintf(ofp, "%lf    %lf    %lf\n", arr[0][x], arr[1][x], arr[2][x]);   
    }


    fclose(ofp);

    printf("Success!\n");

    return 0;
    
}