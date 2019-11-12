/**
 * Solar System
 *
 * This example integrates all planets of the Solar
 * System. The data comes from the NASA HORIZONS system. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
//all data is from nasa horizons, 30th September 2019. Units of AU.

double ss_pos[11][3] = 
{
    {-3.013424684820004E-03 , 7.577400311699641E-03 , 1.522327948994377E-06},
    {-1.744888637875314E-01 ,-4.244459073219732E-01 ,-1.957052267164663E-02},
    {-5.832166811075774E-01 ,-4.227197431546666E-01 , 2.757923523005813E-02},
    { 9.926803557750922E-01 , 1.171071592931596E-01 ,-5.897765002577186E-06},
    { 8.228862161210234E-01 , 6.587996475726317E-01 ,-3.785334833271938E-01},
    {-1.644137782803167E+00 , 2.530657764233049E-01 , 4.541188790127934E-02},
    {-1.699642400881444E-01 ,-5.250743091389841E+00 , 2.557799438953158E-02},
    { 3.337284978938324E+00 ,-9.463581115929795E+00 , 3.169757971579321E-02},
    { 1.643185373589321E+01 , 1.110227970027843E+01 ,-1.716425425492940E-01},
    { 2.917750633262754E+01 ,-6.646446374100315E+00 ,-5.355542983258703E-01},
    { 1.269610583507422E+01 ,-3.141113199578617E+01 ,-3.112775015648088E-01},
};
double ss_vel[11][3] = 
{
    {-8.492322341632166E-06 ,-9.353155515358904E-07 , 2.301541419735291E-07},
    { 2.048459325443986E-02 ,-8.995044409873299E-03 ,-2.614664569098718E-03},
    { 1.189797004491595E-02 ,-1.634153795167905E-02 ,-9.110754924475504E-04},
    {-2.170315646338148E-03 , 1.703520105098581E-02 ,-5.463857304374388E-07},
    {-1.381007577874730E-02 , 5.897333087115376E-03 , 2.754466583425123E-03},
    {-1.556986377304836E-03 ,-1.264431146517457E-02 ,-2.267130777538514E-04},
    { 7.450947804402543E-03 , 1.166544750377484E-04 ,-1.671012875749180E-04},
    { 4.952121119936067E-03 , 1.839073137038874E-03 ,-2.293132844397104E-04},
    {-2.230759206307085E-03 , 3.075630739324861E-03 , 4.027037883636828E-05},
    { 6.759177587143499E-04 , 3.079179855664010E-03 ,-7.882476544965271E-05},
    { 2.988188173523851E-03 , 5.096901737398172E-04 ,-9.289666940024388E-04},
    
};

double ss_mass[11] =
{
    1.988544e30,
    3.302e23,
    48.685e23,
    6.045476731e24,
    1.3e14,
    6.4185e23,
    1898.13e24,
    5.68319e26,
    86.8103e24,
    102.41e24,
    1.4639248e+22,

};


uint32_t hash[11] = {};
double mass[11] = {};
double radius[11] = {};
double xyz[11][3] = {};
double vxvyvz[11][3] = {};
double xyzvxvyvz[11][6] = {};
double com[3] = {};
double lya = 0;
double meg = 0;

void heartbeat(struct reb_simulation* r);
double e_init;
double tmax;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    // Setup constants
    r->dt             = 0.1;                // in days
    tmax            = 7.3e10;            // 200 Myr
    r->G            = 1.4880826e-34;        // in AU^3 / kg / day^2.
    
    //r->integrator        = REB_INTEGRATOR_WHFAST;
    //r->ri_whfast.safe_mode     = 0;        // Turn off safe mode. Need to call reb_integrator_synchronize() before outputs. 
    //r->ri_whfast.corrector     = 11;        // 11th order symplectic corrector
    //r->integrator = REB_INTEGRATOR_JANUS;
    //r->ri_janus.order = 10;
    r->heartbeat        = heartbeat;
    r->exact_finish_time = 1; // Finish exactly at tmax in reb_integrate(). Default is already 1.
    //r->usleep = 10000;
    r->integrator        = REB_INTEGRATOR_IAS15;        // Alternative non-symplectic integrator

    // Initial conditions
    int planets[4] = {0,3,4,5};
    for (int j=0;j<11;j++){
        //int j = planets[i];
        struct reb_particle p = {0};
        p.x  = ss_pos[j][0];         p.y  = ss_pos[j][1];         p.z  = ss_pos[j][2];
        p.vx = ss_vel[j][0];         p.vy = ss_vel[j][1];         p.vz = ss_vel[j][2];
        p.m  = ss_mass[j];
        reb_add(r, p); 
    }



    reb_move_to_com(r);
    e_init = reb_tools_energy(r);
    system("rm -f xyz.txt");
    system("rm -f lyap.txt");
    system("rm -f megno.txt");
    reb_tools_megno_init(r);
    reb_integrate(r, tmax);

}

void heartbeat(struct reb_simulation* r){
    if (reb_output_check(r, 100.)){
        reb_output_timing(r, tmax);
        reb_integrator_synchronize(r);
        //reb_serialize_particle_data(r,hash,mass,radius,xyz,vxvyvz,xyzvxvyvz);
        //com[0] = reb_get_com(r).x;
        //com[1] = reb_get_com(r).y;
        //com[2] = reb_get_com(r).z;
        
        //FILE* f = fopen("xyz.txt","a");
        //double e = reb_tools_energy(r);
        //fprintf(f,"%e %e\n",r->t, fabs((e-e_init)/e_init));
        //fprintf(f,"%e %e %e %e %e %e %e %e %e %e %e %e\n",
        //    com[0],com[1],com[2],
        //    xyz[0][0],xyz[0][1],xyz[0][2],
        //    xyz[3][0],xyz[3][1],xyz[3][2],
        //    xyz[4][0],xyz[4][1],xyz[4][2]);
        //fclose(f);

        lya = reb_tools_calculate_lyapunov(r);
        meg = reb_tools_calculate_megno(r);

        FILE* l = fopen("lyap.txt","a");
        fprintf(l,"%e\n",lya);
        fclose(l);

        FILE* m = fopen("megno.txt","a");
        fprintf(m,"%e\n",meg);
        fclose(m);
        //printf("%f",r->particles[9].x-xyz[9][0]);
        //particles[3] is earth.
        //r->particles[3].x+= -(r->particles[3].x) +ss_pos[3][0];
        //r->particles[3].y+= -(r->particles[3].y) +ss_pos[3][1];
        //r->particles[3].z+= -(r->particles[3].z) +ss_pos[3][2];
    }
}

