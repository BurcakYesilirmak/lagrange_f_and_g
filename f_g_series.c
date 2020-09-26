// This code calculates 2BP using STUMPFF FUNCTION known its parameteres as (Cz, Sz)
// Burcak Yesilirmak


    #include<stdio.h>
    #include<math.h>
    #define mu 398600    // [km^3/s^2] Earth’s gravitational parameter

    // Calculating Stmupff function's Cz and Sz component 
    void Stumpff_Sz(double zi, double *Sz){
         if (zi > 0)
        *Sz = (sqrt(zi) - sin(sqrt(zi)))/pow(sqrt(zi),3);
        else if (zi < 0)
        *Sz = (sinh(sqrt(-zi)) - sqrt(-zi))/pow(sqrt(-zi),3);
        else
        *Sz = 1/6;
    }
    void Stumpff_Cz(double zi, double *Cz){
         if (zi > 0)
        *Cz = (1 - cos(sqrt(zi)))/zi;
        else if (zi < 0)
        *Cz =(cosh(sqrt(-zi)) - 1)/(-zi);
        else
        *Cz = 1/2;
    }
    int main(int argc, char **argv) {
        double dt, h, Vr0, fdot, gdot, c, s, norm_r0, norm_rv0, alfa; 
        double X0, Xi, tol, zi, fX, fdX, eps, dot_rv, norm_v0, Sz, Cz;
        double f, g, R_i, R_j, R_k, r, df, dg, V_i, V_j, V_k;
        // initial position and velocity vectors at the time t0
    double r0[]= {7000, -12124, 0};    // (km)
    double v0[]= {2.667, 4.6210, 0};    // (km/s)
    dt= 3600;  // time interval (s)
    norm_r0 = sqrt(r0[0]* r0[0] + r0[1] * r0[1] + r0[2] * r0[2]);
    norm_v0 = sqrt(v0[0]* v0[0] + v0[1] * v0[1] + v0[2] * v0[2]);
    dot_rv = (r0[0] * v0[0]) + (r0[1] * v0[1]) + (r0[2] * v0[2]);
    Vr0 = dot_rv/norm_r0;
    alfa = 2/norm_r0- norm_v0*norm_v0/mu;

    if  (alfa > 0 )
    printf("The trajectory is an ellipse\n");
    else if (alfa < 0)
    printf("The trajectory is an hyperbola\n");
    else if (alfa == 0)
    printf("The trajectory is an parabola\n");
    printf("alfa= %.6le\n",alfa);

    X0 = sqrt(mu)*fabs(alfa)*dt;  // [km^0.5]  Initial estimate of X0
    Xi = X0;
    tol = 1E-15;  //Tolerance
    while(1) { 
    zi = alfa*Xi*Xi;
      Stumpff_Sz(zi, &Sz);
      Stumpff_Sz(zi, &Cz);
    fX  = norm_r0*Vr0/sqrt(mu)*Xi*Xi*Cz + (1 - alfa*norm_r0)*Xi*Xi*Xi*Sz + norm_r0*Xi -sqrt(mu)*dt;
    fdX = norm_r0*Vr0/sqrt(mu)*Xi*(1 - alfa*Xi*Xi*Sz) + (1 - alfa*norm_r0)*Xi*Xi*Cz + norm_r0;
    eps = fX/fdX;
    Xi = Xi - eps;
    if(fabs(eps) < tol) { 
    break;
    } }
    printf("Universal anomaly X = %4.3f [km^0.5] \n\n",Xi);
   // Lagrange f and g coefficients in terms of the universal anomaly
    f  =  1 - Xi*Xi/norm_r0*Cz;
    printf("f= %4.2f\n",f);
    g  = dt - 1/sqrt(mu)*Xi*Xi*Xi*Sz;
    printf("g= %4.2f\n",g);
    R_i  = f*r0[0] + g*v0[0];
    R_j  = f*r0[1] + g*v0[1];
    R_k  = f*r0[2] + g*v0[2];
    r = sqrt(R_i*R_i+R_j*R_j+R_k*R_k);
    printf("r= %4.2f\n",r);
    // derivatives
    df = sqrt(mu)/(r*norm_r0)*(alfa*Xi*Xi*Xi*Sz - Xi);      
    dg = 1 - Xi*Xi/r*Cz;
    V_i = df*r0[0] + dg*v0[0];
    V_j = df*r0[1] + dg*v0[1];
    V_k = df*r0[2] + dg*v0[2];
    printf("Position after %4.2f [min] R = %4.2f*i + %4.2f*j + %4.2f*k[km] \n",dt/60,R_i,R_j,R_k);
    printf("Velocity after %4.2f [min] V = %4.3f*i + %4.3f*j + %4.2f*k[km/s] \n",dt/60,V_i,V_j,V_k);
    }

