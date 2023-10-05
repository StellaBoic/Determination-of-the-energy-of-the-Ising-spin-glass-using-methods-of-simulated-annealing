#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran1.c" // random number generator

#define max 10   // grid dimension max L max
#define kT0 5.0 // initial temperature
#define dkT 0.02 // temperature step

#define NbSkip 100 // skipped blocks
#define Nk    1000 // number of steps
#define Nb    1100 // number of blocks
#define Nd    1000 // number of steps for the demon algorithm

double energy ( int s[max+2][max+2], int J[2*max+1][max+1])
{
    int i, j, down, right;
    double sum = 0.0;
    for( i=1; i <= max; i++)
        {
            for( j=1; j <= max; j++)
                {
                    // we ensure that even edge elements have a sufficient number of neighbors
                    if(i==max)
                        {
                            right = 1;
                        }
                    else
                        {
                            right = i+1;
                        }
                    if(j==max)
                        {
                            down = 1;
                        }
                    else
                        {
                            down = j+1;
                        }
                    sum += s[i][j]*(interaction( i, i, j, down, J)*s[i][down]+interaction( i, right, j, j, J)*s[right][j]);
                }
        }
    return  - sum;
}

void initial_configuration (int s[max+2][max+2])
{
    long idum=(-2123);
    int i, j;
    for(i=1; i <= max; i++)
        {
            for(j=1; j <= max; j++)
                {
                    //s[i][j] == 1;
                    if(ran1(&idum)>= 0.5)
                        {
                            s[i][j] = 1;
                        }
                    else
                        {
                            s[i][j] = -1;
                        }
                }
        }
}

void setting_configuration (int J[ 2*max+1][max+1])
{
    long idum=(-2123);
    int i, j;
    for(i=1; i <= 2*max; i++)
        {
            for(j=1; j<= max; j++)
                {

                    if(ran1(&idum)>= 0.5)
                        {
                            J[i][j] = 1;
                        }
                    else
                        {
                            J[i][j] = -1;
                        }

                }
        }
}

int interaction (int i1, int i2, int j1, int j2, int J[2*max+1][max+1])
{
    int k,l;

    // boundary conditions
    if(i2 == max+1)
        {
            i2 = 1;
        }
    if(j2==max+1)
        {
            j2 = 1;
        }

    l = i1+i2;
    if(j1>j2)
        {
            k = j1;
        }
    else
        {
            k = j2;
        }
    return J[l][k];
}

double metropolis(int s[max+2][max+2],int J[2*max+1][max+1], double kT, double *reject, FILE *dat1, FILE *dat2, FILE *dat3)
{
    long idum=(-2123);
    int ib, ik, NbEff, i, j;
    double E, dE, sigmaE;
    double SkE, SkE2, SbE, SbE2, aE, aE2;


    SbE = 0.0;
    SbE2 = 0.0;
    *reject = 0.0;

    E = energy(s,J);

    // loops through the blocks
            for( ib = 1; ib <= Nb; ib++)
                {
                    SkE = 0.0;
                    SkE2 = 0.0;

                    // we randomly select an element and change its spin
                    for(ik=1; ik <= Nk;  ik++)
                        {
                            i = (int)(1.+ran1(&idum)*max); // choose an element from 1 to max
                            j = (int)(1.+ran1(&idum)*max);
                            s [i][j] *= -1;
                            // boundary conditions for s
                            if(i==1)
                                {
                                    s[i-1][j] = s[max][j];
                                }
                            if(i==max)
                                {
                                    s[i+1][j] = s[1][j];
                                }
                            if(j==1)
                                {
                                    s[i][j-1] = s[i][max];
                                }
                            if(j==max)
                                {
                                    s[i][j+1] = s[i][1];
                                }
                            // system energy change
                            dE = energy(s,J) - E;

                            // rejection
                            if ( (dE>0) && (exp((-dE)/kT) <= ran1(&idum)) )
                                {
                                    s[i][j] = s[i][j]*(-1);
                                    *reject += 1./Nk/Nb; // percentage of rejected steps
                                    dE = 0.0;
                                }
                            E += dE;


                            if(ib>NbSkip)
                                {
                                    SkE += E;
                                    SkE2 += E*E;
                                }
                        }
                    if (ib>NbSkip)
                        {
                            SbE += SkE/Nk;
                            SbE2 += SkE2/Nk;

                            NbEff = ib-NbSkip;
                            aE = SbE/NbEff;
                            aE2 = SbE2/NbEff;
                            sigmaE = sqrt(aE2 - aE*aE)/sqrt(NbEff);
                            fprintf(dat1, "%d  %lf  %lf  %lf\n", ib, SkE/Nk, aE, sigmaE);
                        }
                }


    // saving the final configuration
    for( i=1; i<= max; i++)
        {
            for( j=1 ; j<=max; j++)
                {
                    if(s[i][j] == 1)
                        {
                            fprintf(dat2,"%lf %d %d \n",kT, i ,j );
                        }
                    if(s[i][j] == -1)
                        {
                            fprintf(dat3,"%lf %d %d \n",kT, i ,j );
                        }
                }
        }
    fprintf(dat2,"\n\n");
    fprintf(dat3,"\n\n");

    return aE/(max*max);
}

double demon(int s[max+2][max+2],int J[2*max+1][max+1], double kT, double *reject, FILE *dat)
{
    long idum=(-2123);
    int i, j, k, counter = 0.0;
    double Es, Ed, dE, sEd= 0.0, sEs = 0.0;

    *reject = 0.0;
    Ed = kT;
    Es = energy(s,J);

    for( k = 0; k < Nd; k++)
        {
            i = (int)(1.+ran1(&idum)*max); // choose an element from 1 to max
            j = (int)(1.+ran1(&idum)*max);
            s [i][j] *= -1;
            // boundary conditions
            if(i==1)
                {
                    s[i-1][j] = s[max][j];
                }
            if(i==max)
                {
                    s[i+1][j] = s[1][j];
                }
            if(j==1)
                {
                    s[i][j-1] = s[i][max];
                }
            if(j==max)
                {
                    s[i][j+1] = s[i][1];
                }
            // system energy change
             dE = energy(s,J) - Es;

           // check to accept the change
            if(Ed >= dE)
                {
                Ed -= dE;
                Es += dE;
                sEd += Ed;
                sEs += Es;
                counter++;
                }
            else
                {
                s [i][j] *= -1;
                *reject += 1./Nd;
                }
        }

    sEd /= counter;
    sEs /= counter;

    return sEs/(max*max);
}

double DkT( int it)
{
    return exp( - 5.0/it);
}

int main()
{
    long idum=(-2123);
    int s[max+2][max+2], J[2*max+1][max+1], i, j, k = 0.0;
    double kT = kT0, E, E0 = 1000.0, T0, reject, reject0 = 0.0;
    int it, ib, ik, NbEff;

    // files
    FILE *dat1, *dat2, *dat3, *dat4;
    dat1 = fopen("Demon.Energy.lin.dat","w");
    dat2 = fopen("Demon.E-T-r.lin.dat","w");
    dat3 = fopen("Demon.Spin_up.lin.dat","w");
    dat4 = fopen("Demon.Spin_down.lin.dat","w");

    FILE *dat5, *dat6, *dat7, *dat8;
    dat5 = fopen("Demon.Energy.dat","w");
    dat6 = fopen("Demon.E-T-r.dat","w");
    dat7 = fopen("Demon.Spin_up.dat","w");
    dat8 = fopen("Demon.Spin_down.dat","w");



    // we set the matrix of spins and take the random matrix of interactions Jij
    initial_configuration(s);
    setting_configuration(J);

    // LINEAR SCHEDULE OF COOLING
    kT = kT0;
    E0 = 1000.0;
    for(it=0; kT > 0 ; it++)
        {
            kT = kT0 - it*dkT;

            // E = metropolis(s,J,kT,&reject,dat1,dat3,dat4);
            E = demon(s,J,kT,&reject,dat1);

            fprintf(dat2,"%lf\t%lf\t%lf\n", kT, reject, E);
            printf("kT = %lf  E= %lf \n",kT, E);

            // we are looking for minimum energy and corresponding temperature
            if(E < E0)
                {
                    E0 = E;
                    T0 = kT;
                }
        }
    printf("E0/N = %lf Tss=%lf\n",E0,T0);


    // EXPONENTIAL COOLING SCHEDULE

    initial_configuration(s);

    kT = kT0;
    E0 = 1000.0;
    for(it=0; kT > 0 ; it++)
        {
            kT = kT0 - DkT(it)*kT0 - 0.1; // exponential cooling schedule

            // E = metropolis(s,J,kT,&reject,dat5,dat7,dat8);
            E = demon(s,J,kT,&reject,dat5);

            fprintf(dat6,"%lf\t%lf\t%lf\n", kT, reject, E);
            printf("kT = %lf  E= %lf \n",kT, E);

            // we are looking for minimum energy and corresponding temperature
            if(E < E0 && kT > 0)
                {
                    E0 = E;
                    T0 = kT;
                }
        }
    printf("E0/N = %lf Tss=%lf\n",E0,T0);


    fclose(dat1);
    fclose(dat2);
    fclose(dat3);
    fclose(dat4);
    fclose(dat5);
    fclose(dat6);
    fclose(dat7);
    fclose(dat8);
    return 0;
}
