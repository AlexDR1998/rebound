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
double semi_major = 0;
double ml_earth =0;
double ml_cruithne=0;
double ve_cr_dist_sq = 0;
double arg_per_cruithne=0;
double inc_cruithne=0;
double l_asc_cruithne=0;
double tot_ang_mom[3]={};
double eccentricity = 0;
double mean_anom=0;
int reversible_test=1;

void heartbeat(struct reb_simulation* r);
void heartbeat2(struct reb_simulation* r);
double e_init;
double tmax;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    // Setup constants
    r->dt             = 1;                // in days
    tmax            = 2500000;//7.3e10;            // 200 Myr
    //tmax              = 1000000;
    r->G            = 1.4880826e-34;        // in AU^3 / kg / day^2.
    
    //r->integrator        = REB_INTEGRATOR_WHFAST;
    //r->ri_whfast.safe_mode     = 0;        // Turn off safe mode. Need to call reb_integrator_synchronize() before outputs. 
    //r->ri_whfast.corrector     = 17;        // 17th order symplectic corrector
    r->integrator = REB_INTEGRATOR_JANUS;
    r->ri_janus.order = 8;
    r->ri_janus.scale_pos = 1e-17;
    r->ri_janus.scale_vel = 1e-17*3.14159265359;
    r->ri_janus.recalculate_integer_coordinates_this_timestep=1;

    r->heartbeat        = heartbeat;
    r->exact_finish_time = 1; // Finish exactly at tmax in reb_integrate(). Default is already 1.
    //r->usleep = 10000;
    //r->integrator        = REB_INTEGRATOR_IAS15;        // Alternative non-symplectic integrator

    // Initial conditions
    //double dir=-1;
    for (int i=0;i<11;i++){
        struct reb_particle p = {0};
        p.x  = ss_pos[i][0];         p.y  = ss_pos[i][1];         p.z  = ss_pos[i][2];
        p.vx = ss_vel[i][0];         p.vy = ss_vel[i][1];         p.vz = ss_vel[i][2];
        p.m  = ss_mass[i];
        reb_add(r, p); 
    }



    reb_move_to_com(r);
    e_init = reb_tools_energy(r);
    system("rm -f xyz.txt");
    system("rm -f cruithne_semi_major_axis_janus_rev_ord8.txt");
    system("rm -f longitudes_janus_rev_ord8.txt");
    system("rm -f energies_janus_rev_ord8.txt");
    system("rm -f ven_dist.txt");
    system("rm -f inclinations_janus_rev_ord8.txt");
    system("rm -f longitude_ascending_janus_rev_ord8.txt");
    system("rm -f argument_peri_janus_rev_ord8.txt");
    system("rm -f angular_momentum_janus_rev_ord8.txt");
    system("rm -f eccentricity_janus_rev_ord8.txt");
    system("rm -f mean_anom_janus_rev_ord8.txt");
    reb_integrate(r, tmax);


    //Pause and rerun new simulation with speeds reversed
    struct reb_simulation* r2 = reb_create_simulation();

    r2->dt=1;
    r2->G=1.4880826e-34;
    r2->heartbeat        = heartbeat2;
    r2->exact_finish_time = 1;


    r2->integrator = REB_INTEGRATOR_JANUS;
    r2->ri_janus.order = 8;
    r2->ri_janus.scale_pos = 1e-17;
    r2->ri_janus.scale_vel = 1e-17*3.14159265359;
    r2->ri_janus.recalculate_integer_coordinates_this_timestep=1;
    //r2->integrator        = REB_INTEGRATOR_IAS15;        // Alternative non-symplectic integrator
    //r2->integrator        = REB_INTEGRATOR_WHFAST;
    //r2->ri_whfast.safe_mode     = 0;        // Turn off safe mode. Need to call reb_integrator_synchronize() before outputs. 
    //r2->ri_whfast.corrector     = 17;        // 17th order symplectic corrector
    for (int i=0;i<11;i++){
        //for (int j=0;j<3;j++){
        //    vxvyvz[i][j] =-1*vxvyvz[i][j];
        //    xyzvxvyvz[i][j+3] = -1*xyzvxvyvz[i][j+3];
        //}
        struct reb_particle p = {0};
        p.x  = xyz[i][0];         p.y  = xyz[i][1];         p.z  = xyz[i][2];
        p.vx = -vxvyvz[i][0];     p.vy = -vxvyvz[i][1];     p.vz = -vxvyvz[i][2];
        p.m  = ss_mass[i];
        reb_add(r2, p); 

    }

    //reb_set_serialized_particle_data(r2,hash,mass,radius,xyz,vxvyvz,xyzvxvyvz);
    reb_move_to_com(r2);
    reb_integrate(r2, tmax);

}

void heartbeat(struct reb_simulation* r){
    


    if (reb_output_check(r, 10.)){
        reb_output_timing(r, tmax);
        reb_integrator_synchronize(r);


        //---Get data out of simulation        
        reb_serialize_particle_data(r,hash,mass,radius,xyz,vxvyvz,xyzvxvyvz);
        semi_major  = reb_tools_particle_to_orbit(r->G,r->particles[4],r->particles[0]).a;
        ml_earth    = reb_tools_particle_to_orbit(r->G,r->particles[3],r->particles[0]).l;
        ml_cruithne = reb_tools_particle_to_orbit(r->G,r->particles[4],r->particles[0]).l;
        inc_cruithne= reb_tools_particle_to_orbit(r->G,r->particles[4],r->particles[0]).inc;
        l_asc_cruithne=reb_tools_particle_to_orbit(r->G,r->particles[4],r->particles[0]).Omega;
        arg_per_cruithne=reb_tools_particle_to_orbit(r->G,r->particles[4],r->particles[0]).omega;
        eccentricity=reb_tools_particle_to_orbit(r->G,r->particles[4],r->particles[0]).e;
        mean_anom=reb_tools_particle_to_orbit(r->G,r->particles[4],r->particles[0]).M;

        //com[0] = reb_get_com(r).x;
        //com[1] = reb_get_com(r).y;
        //com[2] = reb_get_com(r).z;
        //tot_ang_mom=0;
        //for(int i=1;i<11;i++){
        //    tot_ang_mom+=reb_tools_particle_to_orbit(r->G,r->particles[i],r->particles[0]).h;
        //}



        //---Write data to txt files
        //FILE* f = fopen("xyz.txt","a");
        //fprintf(f,"%e %e %e %e %e %e %e %e %e %e %e %e\n",
        //    com[0],com[1],com[2],
        //    xyz[0][0],xyz[0][1],xyz[0][2],
        //    xyz[3][0],xyz[3][1],xyz[3][2],
        //    xyz[4][0],xyz[4][1],xyz[4][2]);
        //fclose(f);
        //Distance between venus and Cruithne
        //ve_cr_dist_sq=0;
        //for(int i=0;i<3;i++){
        //    ve_cr_dist_sq+=(xyz[2][i]-xyz[4][i])*(xyz[2][i]-xyz[4][i]);
        //}
        //FILE* ven_dist = fopen("ven_dist.txt","a");
        //fprintf(ven_dist, "%e\n", sqrt(ve_cr_dist_sq));
        //fclose(ven_dist);

        FILE* sm = fopen("cruithne_semi_major_axis_janus_rev_ord8.txt","a");
        fprintf(sm,"%e\n",semi_major);
        fclose(sm);

        FILE* inc = fopen("inclinations_janus_rev_ord8.txt","a");
        fprintf(inc, "%e\n",inc_cruithne);
        fclose(inc);

        FILE* asc = fopen("longitude_ascending_janus_rev_ord8.txt","a");
        fprintf(asc, "%e\n",l_asc_cruithne);
        fclose(asc);

        FILE* arg = fopen("argument_peri_janus_rev_ord8.txt","a");
        fprintf(arg, "%e\n",arg_per_cruithne);
        fclose(arg);

        //FILE* ang = fopen("angular_momentum_janus_rev_ord8.txt","a");
        //fprintf(ang, "%e\n",);

        FILE* ml = fopen("longitudes_janus_rev_ord8.txt","a");
        fprintf(ml,"%e %e\n", ml_earth, ml_cruithne);
        fclose(ml);

        FILE* man=fopen("mean_anom_janus_rev_ord8.txt","a");
        fprintf(man, "%e\n",mean_anom);
        fclose(man);

        FILE* ecc=fopen("eccentricity_janus_rev_ord8.txt","a");
        fprintf(ecc, "%e\n",eccentricity);
        fclose(ecc);
    }



}



void heartbeat2(struct reb_simulation* r2){
    


    if (reb_output_check(r2, 10.)){
        reb_output_timing(r2, tmax);
        reb_integrator_synchronize(r2);


        //---Get data out of simulation        
        reb_serialize_particle_data(r2,hash,mass,radius,xyz,vxvyvz,xyzvxvyvz);
        semi_major  = reb_tools_particle_to_orbit(r2->G,r2->particles[4],r2->particles[0]).a;
        ml_earth    = reb_tools_particle_to_orbit(r2->G,r2->particles[3],r2->particles[0]).l;
        ml_cruithne = reb_tools_particle_to_orbit(r2->G,r2->particles[4],r2->particles[0]).l;
        inc_cruithne= reb_tools_particle_to_orbit(r2->G,r2->particles[4],r2->particles[0]).inc;
        l_asc_cruithne=reb_tools_particle_to_orbit(r2->G,r2->particles[4],r2->particles[0]).Omega;
        arg_per_cruithne=reb_tools_particle_to_orbit(r2->G,r2->particles[4],r2->particles[0]).omega;
        eccentricity=reb_tools_particle_to_orbit(r2->G,r2->particles[4],r2->particles[0]).e;
        mean_anom=reb_tools_particle_to_orbit(r2->G,r2->particles[4],r2->particles[0]).M;

        //com[0] = reb_get_com(r).x;
        //com[1] = reb_get_com(r).y;
        //com[2] = reb_get_com(r).z;
        //tot_ang_mom=0;
        //for(int i=1;i<11;i++){
        //    tot_ang_mom+=reb_tools_particle_to_orbit(r->G,r->particles[i],r->particles[0]).h;
        //}



        //---Write data to txt files
        //FILE* f = fopen("xyz.txt","a");
        //fprintf(f,"%e %e %e %e %e %e %e %e %e %e %e %e\n",
        //    com[0],com[1],com[2],
        //    xyz[0][0],xyz[0][1],xyz[0][2],
        //    xyz[3][0],xyz[3][1],xyz[3][2],
        //    xyz[4][0],xyz[4][1],xyz[4][2]);
        //fclose(f);
        //Distance between venus and Cruithne
        //ve_cr_dist_sq=0;
        //for(int i=0;i<3;i++){
        //    ve_cr_dist_sq+=(xyz[2][i]-xyz[4][i])*(xyz[2][i]-xyz[4][i]);
        //}
        //FILE* ven_dist = fopen("ven_dist.txt","a");
        //fprintf(ven_dist, "%e\n", sqrt(ve_cr_dist_sq));
        //fclose(ven_dist);

        FILE* sm = fopen("cruithne_semi_major_axis_janus_rev_ord8.txt","a");
        fprintf(sm,"%e\n",semi_major);
        fclose(sm);

        FILE* inc = fopen("inclinations_janus_rev_ord8.txt","a");
        fprintf(inc, "%e\n",inc_cruithne);
        fclose(inc);

        FILE* asc = fopen("longitude_ascending_janus_rev_ord8.txt","a");
        fprintf(asc, "%e\n",l_asc_cruithne);
        fclose(asc);

        FILE* arg = fopen("argument_peri_janus_rev_ord8.txt","a");
        fprintf(arg, "%e\n",arg_per_cruithne);
        fclose(arg);

        //FILE* ang = fopen("angular_momentum_janus_rev_ord8.txt","a");
        //fprintf(ang, "%e\n",);

        FILE* ml = fopen("longitudes_janus_rev_ord8.txt","a");
        fprintf(ml,"%e %e\n", ml_earth, ml_cruithne);
        fclose(ml);

        FILE* man=fopen("mean_anom_janus_rev_ord8.txt","a");
        fprintf(man, "%e\n",mean_anom);
        fclose(man);

        FILE* ecc=fopen("eccentricity_janus_rev_ord8.txt","a");
        fprintf(ecc, "%e\n",eccentricity);
        fclose(ecc);
    }


}
