/* ========================================================================= *
 * vehicle_14dof_mex.c
 * Continuous C-MEX S-Function for 14-DOF Vehicle Dynamics & Non-linear Tire
 * Integrated with Steer-by-Wire (SbW) Front Wheel Steering Dynamics
 * Equations (19) to (43) of Nugroho et al., IEEE Trans. Reliability
 * ========================================================================= */

#define S_FUNCTION_NAME  vehicle_14dof_mex
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <math.h>

#define NUM_STATES 28
#define NUM_INPUTS 3   /* u[0] = delta_f (rad), u[1] = T_drive (Nm), u[2] = T_brake (Nm) */
#define NUM_OUTPUTS 8  /* u, v, r, phi, p, theta, q, Tz */

/* Vehicle Physical Parameters */
#define MASS_TOTAL    1500.0   /* m (kg) */
#define MASS_SPRUNG   1350.0   /* m_s (kg) */
#define MASS_UNSPRUNG 37.5     /* m_u (kg) per wheel */
#define I_XX          550.0    /* kg*m^2 */
#define I_YY          1800.0   /* kg*m^2 */
#define I_ZZ          2200.0   /* kg*m^2 */
#define I_WHEEL       1.2      /* kg*m^2 */
#define TRACK_F       1.55     /* t_f (m) */
#define TRACK_R       1.53     /* t_r (m) */
#define DIST_A        1.15     /* a (m) */
#define DIST_B        1.45     /* b (m) */
#define H_CG          0.52     /* h_cg (m) */
#define R_W           0.31     /* Wheel radius (m) */
#define K_SUSP        35000.0  /* Suspension stiffness (N/m) */
#define B_SUSP        2500.0   /* Suspension damping (N*s/m) */
#define GRAVITY       9.81     /* g (m/s^2) */
#define TRAIL_TP      0.035    /* Pneumatic trail (m) */
#define TRAIL_TM      0.015    /* Mechanical castor trail (m) */

/* Pacejka Magic Formula Coefficients */
#define PAC_B  10.0
#define PAC_C  1.30
#define PAC_D  1.00
#define PAC_E  -0.8

static double compute_pacejka_Fy(double alpha, double Fz) {
    double B_al = PAC_B * alpha;
    double phi_val = B_al - PAC_E * (B_al - atan(B_al));
    return Fz * PAC_D * sin(PAC_C * atan(phi_val));
}

static void mdlInitializeSizes(SimStruct *S) {
    ssSetNumSFcnParams(S, 0);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) return;

    ssSetNumContStates(S, NUM_STATES);
    ssSetNumDiscStates(S, 0);

    if (!ssSetNumInputPorts(S, 1)) return;
    ssSetInputPortWidth(S, 0, NUM_INPUTS);
    ssSetInputPortDirectFeedThrough(S, 0, 1);

    if (!ssSetNumOutputPorts(S, 1)) return;
    ssSetOutputPortWidth(S, 0, NUM_OUTPUTS);

    ssSetNumSampleTimes(S, 1);
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}

static void mdlInitializeSampleTimes(SimStruct *S) {
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_INITIALIZE_CONDITIONS
static void mdlInitializeConditions(SimStruct *S) {
    real_T *x0 = ssGetContStates(S);
    int_T i;
    for (i = 0; i < NUM_STATES; i++) {
        x0[i] = 0.0;
    }
    x0[1] = 20.0; /* Initial longitudinal velocity u = 20 m/s (72 km/h) */
    x0[24] = 20.0 / R_W; /* fl spin */
    x0[25] = 20.0 / R_W; /* fr spin */
    x0[26] = 20.0 / R_W; /* rl spin */
    x0[27] = 20.0 / R_W; /* rr spin */
}

static void mdlOutputs(SimStruct *S, int_T tid) {
    real_T *y = ssGetOutputPortRealSignal(S, 0);
    real_T *x = ssGetContStates(S);
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S, 0);
    double delta_f = *uPtrs[0];

    /* States */
    double u_vel = x[1];
    double v_vel = x[3];
    double phi   = x[6];
    double p     = x[7];
    double theta = x[8];
    double q     = x[9];
    double r     = x[11];

    /* Slip angles contact patch (Pers. 38 - 39) */
    double denom_fl = u_vel - 0.5 * TRACK_F * r;
    double denom_fr = u_vel + 0.5 * TRACK_F * r;
    if (fabs(denom_fl) < 0.1) denom_fl = 0.1;
    if (fabs(denom_fr) < 0.1) denom_fr = 0.1;

    double alpha_fl = delta_f - atan((v_vel + DIST_A * r) / denom_fl);
    double alpha_fr = delta_f - atan((v_vel + DIST_A * r) / denom_fr);

    /* Normal load static & dynamic transfer (Pers. 41) */
    double Fz_static_f = MASS_TOTAL * GRAVITY * (DIST_B / (DIST_A + DIST_B)) * 0.5;
    double Fy_fl = compute_pacejka_Fy(alpha_fl, Fz_static_f);
    double Fy_fr = compute_pacejka_Fy(alpha_fr, Fz_static_f);

    /* Cumulative self-aligning torque Tz (Pers. 42 - 43) */
    double Tz = (Fy_fl + Fy_fr) * (TRAIL_TP + TRAIL_TM);

    y[0] = u_vel;
    y[1] = v_vel;
    y[2] = r;
    y[3] = phi;
    y[4] = p;
    y[5] = theta;
    y[6] = q;
    y[7] = Tz;
}

#define MDL_DERIVATIVES
static void mdlDerivatives(SimStruct *S) {
    real_T *dx = ssGetdX(S);
    real_T *x  = ssGetContStates(S);
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S, 0);

    double delta_f = *uPtrs[0];
    double T_drive = *uPtrs[1];
    double T_brake = *uPtrs[2];

    /* Unpack Chassis States */
    double u_vel = x[1];
    double v_vel = x[3];
    double ws    = x[5];
    double phi   = x[6];
    double p     = x[7];
    double theta = x[8];
    double q     = x[9];
    double psi   = x[10];
    double r     = x[11];

    /* Unsprung displacements & velocities */
    double zu_fl = x[12], zu_fr = x[13], zu_rl = x[14], zu_rr = x[15];
    double dzu_fl = x[16], dzu_fr = x[17], dzu_rl = x[18], dzu_rr = x[19];

    /* Wheel speeds */
    double om_fl = x[24], om_fr = x[25], om_rl = x[26], om_rr = x[27];

    /* Corner displacements & velocities (Pers. 35 - 36) */
    double zs_fl = x[4] + 0.5 * TRACK_F * phi - DIST_A * theta;
    double zs_fr = x[4] - 0.5 * TRACK_F * phi - DIST_A * theta;
    double zs_rl = x[4] + 0.5 * TRACK_R * phi + DIST_B * theta;
    double zs_rr = x[4] - 0.5 * TRACK_R * phi + DIST_B * theta;

    double dzs_fl = ws + 0.5 * TRACK_F * p - DIST_A * q;
    double dzs_fr = ws - 0.5 * TRACK_F * p - DIST_A * q;
    double dzs_rl = ws + 0.5 * TRACK_R * p + DIST_B * q;
    double dzs_rr = ws - 0.5 * TRACK_R * p + DIST_B * q;

    /* Suspension forces (Pers. 34) */
    double Fsusp_fl = K_SUSP * (zs_fl - zu_fl) + B_SUSP * (dzs_fl - dzu_fl);
    double Fsusp_fr = K_SUSP * (zs_fr - zu_fr) + B_SUSP * (dzs_fr - dzu_fr);
    double Fsusp_rl = K_SUSP * (zs_rl - zu_rl) + B_SUSP * (dzs_rl - dzu_rl);
    double Fsusp_rr = K_SUSP * (zs_rr - zu_rr) + B_SUSP * (dzs_rr - dzu_rr);

    /* Contact Normal Loads Fz */
    double Fz_fl = Fsusp_fl + MASS_UNSPRUNG * GRAVITY;
    double Fz_fr = Fsusp_fr + MASS_UNSPRUNG * GRAVITY;
    double Fz_rl = Fsusp_rl + MASS_UNSPRUNG * GRAVITY;
    double Fz_rr = Fsusp_rr + MASS_UNSPRUNG * GRAVITY;
    if (Fz_fl < 0.0) Fz_fl = 0.0;
    if (Fz_fr < 0.0) Fz_fr = 0.0;
    if (Fz_rl < 0.0) Fz_rl = 0.0;
    if (Fz_rr < 0.0) Fz_rr = 0.0;

    /* Slip angles contact patch (Pers. 38 - 39) */
    double d_fl = u_vel - 0.5 * TRACK_F * r; if (fabs(d_fl) < 0.1) d_fl = 0.1;
    double d_fr = u_vel + 0.5 * TRACK_F * r; if (fabs(d_fr) < 0.1) d_fr = 0.1;
    double d_rl = u_vel - 0.5 * TRACK_R * r; if (fabs(d_rl) < 0.1) d_rl = 0.1;
    double d_rr = u_vel + 0.5 * TRACK_R * r; if (fabs(d_rr) < 0.1) d_rr = 0.1;

    double alpha_fl = delta_f - atan((v_vel + DIST_A * r) / d_fl);
    double alpha_fr = delta_f - atan((v_vel + DIST_A * r) / d_fr);
    double alpha_rl = -atan((v_vel - DIST_B * r) / d_rl);
    double alpha_rr = -atan((v_vel - DIST_B * r) / d_rr);

    /* Lateral forces */
    double Fy_fl = compute_pacejka_Fy(alpha_fl, Fz_fl);
    double Fy_fr = compute_pacejka_Fy(alpha_fr, Fz_fr);
    double Fy_rl = compute_pacejka_Fy(alpha_rl, Fz_rl);
    double Fy_rr = compute_pacejka_Fy(alpha_rr, Fz_rr);

    /* Longitudinal tire forces approx */
    double Fx_fl = 0.0, Fx_fr = 0.0, Fx_rl = 0.0, Fx_rr = 0.0;
    if (u_vel > 0.5) {
        Fx_fl = (T_drive * 0.25 - T_brake * 0.25) / R_W;
        Fx_fr = Fx_fl; Fx_rl = Fx_fl; Fx_rr = Fx_fl;
    }

    /* Transform to local body frame */
    double Fx_body_fl = Fx_fl * cos(delta_f) - Fy_fl * sin(delta_f);
    double Fy_body_fl = Fx_fl * sin(delta_f) + Fy_fl * cos(delta_f);
    double Fx_body_fr = Fx_fr * cos(delta_f) - Fy_fr * sin(delta_f);
    double Fy_body_fr = Fx_fr * sin(delta_f) + Fy_fr * cos(delta_f);

    double sum_Fx = Fx_body_fl + Fx_body_fr + Fx_rl + Fx_rr;
    double sum_Fy = Fy_body_fl + Fy_body_fr + Fy_rl + Fy_rr;

    /* Total Yaw Moment Myaw (Pers. 32) */
    double Myaw = DIST_A * (Fy_body_fl + Fy_body_fr) - DIST_B * (Fy_rl + Fy_rr)
                + 0.5 * TRACK_F * (Fx_body_fr - Fx_body_fl)
                + 0.5 * TRACK_R * (Fx_rr - Fx_rl);

    /* Chassis State Derivatives (Pers. 20 - 31) */
    dx[0] = u_vel * cos(theta) * cos(psi); /* X dot */
    dx[1] = (1.0 / MASS_TOTAL) * sum_Fx + v_vel * r - ws * q - GRAVITY * sin(theta); /* u dot */
    dx[2] = u_vel * cos(theta) * sin(psi); /* Y dot */
    dx[3] = (1.0 / MASS_TOTAL) * sum_Fy - u_vel * r + ws * p + GRAVITY * sin(phi) * cos(theta); /* v dot */
    dx[4] = ws; /* zs dot */
    dx[5] = (1.0 / MASS_SPRUNG) * (-Fsusp_fl - Fsusp_fr - Fsusp_rl - Fsusp_rr) - GRAVITY; /* ws dot */

    dx[6] = p + q * sin(phi) * tan(theta) + r * cos(phi) * tan(theta); /* phi dot */
    dx[7] = (1.0 / I_XX) * (0.5 * TRACK_F * (Fsusp_fl - Fsusp_fr) + 0.5 * TRACK_R * (Fsusp_rl - Fsusp_rr)
            + MASS_SPRUNG * GRAVITY * H_CG * sin(phi)) + ((I_YY - I_ZZ) / I_XX) * q * r; /* p dot */

    dx[8] = q * cos(phi) - r * sin(phi); /* theta dot */
    dx[9] = (1.0 / I_YY) * (-DIST_A * (Fsusp_fl + Fsusp_fr) + DIST_B * (Fsusp_rl + Fsusp_rr)
            - MASS_SPRUNG * GRAVITY * H_CG * sin(theta)) + ((I_ZZ - I_XX) / I_YY) * r * p; /* q dot */

    dx[10] = q * (sin(phi) / cos(theta)) + r * (cos(phi) / cos(theta)); /* psi dot */
    dx[11] = (Myaw / I_ZZ) + ((I_XX - I_YY) / I_ZZ) * p * q; /* r dot */

    /* Unsprung masses vertical dynamics (Pers. 33) */
    dx[12] = dzu_fl; dx[16] = (1.0 / MASS_UNSPRUNG) * (Fz_fl - Fsusp_fl - MASS_UNSPRUNG * GRAVITY);
    dx[13] = dzu_fr; dx[17] = (1.0 / MASS_UNSPRUNG) * (Fz_fr - Fsusp_fr - MASS_UNSPRUNG * GRAVITY);
    dx[14] = dzu_rl; dx[18] = (1.0 / MASS_UNSPRUNG) * (Fz_rl - Fsusp_rl - MASS_UNSPRUNG * GRAVITY);
    dx[15] = dzu_rr; dx[19] = (1.0 / MASS_UNSPRUNG) * (Fz_rr - Fsusp_rr - MASS_UNSPRUNG * GRAVITY);

    /* Wheel displacements & velocities */
    dx[20] = om_fl; dx[24] = (1.0 / I_WHEEL) * (T_drive * 0.25 - T_brake * 0.25 - R_W * Fx_fl);
    dx[21] = om_fr; dx[25] = (1.0 / I_WHEEL) * (T_drive * 0.25 - T_brake * 0.25 - R_W * Fx_fr);
    dx[22] = om_rl; dx[26] = (1.0 / I_WHEEL) * (T_drive * 0.25 - T_brake * 0.25 - R_W * Fx_rl);
    dx[23] = om_rr; dx[27] = (1.0 / I_WHEEL) * (T_drive * 0.25 - T_brake * 0.25 - R_W * Fx_rr);
}

static void mdlTerminate(SimStruct *S) { }

#ifdef MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif