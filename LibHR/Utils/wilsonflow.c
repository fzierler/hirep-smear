/* arXiv:1006.4518 [hep-lat] */

#include "global.h"
#include "geometry.h"
#include "suN.h"
#include "suN_repr_func.h"
#include "memory.h"
#include "global.h"
#include "logger.h"
#include "update.h"
#include "utils.h"
#include "communications.h"
#include "wilsonflow.h"
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

static suNg_field *ws_gf = NULL;
static suNg_field *ws_gf_tmp = NULL;
static suNg_field *Vprime = NULL;
static suNg_field *u_gauge_backup = NULL;
static double *wf_plaq_weight = NULL;

void WF_initialize()
{
#ifdef BC_XYZ_TWISTED
  error(0 == 0, 0, "WF_initialize", "WF has not yet been setup to work with BC_XYZ_TWISTED enabled");
#endif

  if (ws_gf == NULL)
  {
    ws_gf = alloc_gfield(&glattice);
    ws_gf_tmp = alloc_gfield(&glattice);
    Vprime = alloc_gfield(&glattice);
    u_gauge_backup = alloc_gfield(&glattice);

#ifdef BC_T_OPEN
    int ix, iy, iz, index, mu, nu;
    wf_plaq_weight = malloc(sizeof(double) * glattice.gsize_gauge * 16);
    init_plaq_open_BCs(wf_plaq_weight, NULL, 1.0, 1.0);

    if (COORD[0] == 0)
    {
      for (ix = 0; ix < X_EXT; ++ix)
        for (iy = 0; iy < Y_EXT; ++iy)
          for (iz = 0; iz < Z_EXT; ++iz)
          {
            index = ipt_ext(T_BORDER, ix, iy, iz);
            if (index != -1)
            {

              for (mu = 0; mu < 3; mu++)
                for (nu = mu + 1; nu < 4; nu++)
                {
                  wf_plaq_weight[index * 16 + mu * 4 + nu] *= 0.5;
                  wf_plaq_weight[index * 16 + nu * 4 + mu] *= 0.5;
                }
            }
          }
      if (T_BORDER > 0)
      {
        for (ix = 0; ix < X_EXT; ++ix)
          for (iy = 0; iy < Y_EXT; ++iy)
            for (iz = 0; iz < Z_EXT; ++iz)
            {
              index = ipt_ext(T_BORDER - 1, ix, iy, iz);
              if (index != -1)
              {

                for (mu = 1; mu < 3; mu++)
                  for (nu = mu + 1; nu < 4; nu++)
                  {
                    wf_plaq_weight[index * 16 + mu * 4 + nu] *= 0.5;
                    wf_plaq_weight[index * 16 + nu * 4 + mu] *= 0.5;
                  }
              }
            }
      }
    }

    if (COORD[0] == NP_T - 1)
    {
      if (T_BORDER > 0)
      {
        for (ix = 0; ix < X_EXT; ++ix)
          for (iy = 0; iy < Y_EXT; ++iy)
            for (iz = 0; iz < Z_EXT; ++iz)
            {
              index = ipt_ext(T + T_BORDER, ix, iy, iz);
              if (index != -1)
              {
                for (mu = 0; mu < 3; mu++)
                  for (nu = mu + 1; nu < 4; nu++)
                  {
                    wf_plaq_weight[index * 16 + mu * 4 + nu] *= 0.5;
                    wf_plaq_weight[index * 16 + nu * 4 + mu] *= 0.5;
                  }
              }
            }
      }

      for (ix = 0; ix < X_EXT; ++ix)
        for (iy = 0; iy < Y_EXT; ++iy)
          for (iz = 0; iz < Z_EXT; ++iz)
          {
            index = ipt_ext(T + T_BORDER - 2, ix, iy, iz);
            if (index != -1)
            {
              mu = 0;
              for (nu = mu + 1; nu < 4; nu++)
              {
                wf_plaq_weight[index * 16 + mu * 4 + nu] *= 0.5;
                wf_plaq_weight[index * 16 + nu * 4 + mu] *= 0.5;
              }
            }

            index = ipt_ext(T + T_BORDER - 1, ix, iy, iz);
            if (index != -1)
            {
              for (mu = 1; mu < 3; mu++)
                for (nu = mu + 1; nu < 4; nu++)
                {
                  wf_plaq_weight[index * 16 + mu * 4 + nu] *= 0.5;
                  wf_plaq_weight[index * 16 + nu * 4 + mu] *= 0.5;
                }
            }
          }
    }
#endif
#if defined(PLAQ_WEIGHTS) && !defined(BC_T_OPEN)
    wf_plaq_weight = plaq_weight;
#endif
  }
}

void WF_set_bare_anisotropy(double *wf_chi)
{

#ifdef PURE_GAUGE_ANISOTROPY
  WF_initialize();
  int ix, iy, iz, it, index, mu, nu;
  if (wf_plaq_weight == NULL)
  {
    wf_plaq_weight = malloc(sizeof(double) * glattice.gsize_gauge * 16);
    for (index = 0; index < glattice.gsize_gauge * 16; index++)
      wf_plaq_weight[index] = 1.0;
  }

  for (it = 0; it < T_EXT; ++it)
    for (ix = 0; ix < X_EXT; ++ix)
      for (iy = 0; iy < Y_EXT; ++iy)
        for (iz = 0; iz < Z_EXT; ++iz)
        {
          index = ipt_ext(it, ix, iy, iz);
          for (mu = 1; mu < 3; mu++)
            for (nu = mu + 1; nu < 4; nu++)
            {
              if (index != -1)
              {
                wf_plaq_weight[index * 16 + mu * 4 + nu] *= *wf_chi;
                wf_plaq_weight[index * 16 + nu * 4 + mu] *= *wf_chi;
              }
            }
        }
#else
  error(0 == 0, 0, "WF_set_bare_anisotropy", "In order to use anisotropic lattice you must compile with PURE_GAUGE_ANISOTROPY enabled");
#endif
}

void WF_free()
{
  if (ws_gf != NULL)
  {
    free_gfield(ws_gf);
    free_gfield(ws_gf_tmp);
    free_gfield(Vprime);
    free_gfield(u_gauge_backup);
    if (wf_plaq_weight != NULL)
    {
#ifdef PLAQ_WEIGHTS
      if (wf_plaq_weight != plaq_weight)
        free(wf_plaq_weight);
#else
      free(wf_plaq_weight);
#endif
    }
  }
}

/*
 * d/dt V = Z(V) V
 * S_W = 1/g^2 \sum_{p oriented} Re tr ( 1 - V(p) )
 * d_L f(V) = T^a d/ds f(e^{sT^a} V)
 * T_a^dag = -T_a      tr T_a T_b = -delta_{ab}/2
 * Z(V) = -g^2 d_L S_W = \sum_{p oriented} d_L Re tr ( V(p) )
 * Z(V) = d_L Re ( tr V s + tr V^dag s^dag ) =
 *      = d_L ( tr V s + tr V^dag s^dag ) =
 *      = T^a tr T^a ( V s - s^dag V^dag )
 *      = - 1/2 ( V s - s^dag V^dag - 1/N tr (V s - s^dag V^dag) )
 */
static void Zeta(suNg_field *Z, const suNg_field *U, const double alpha)
{

  error(Z->type != &glattice, 1, "wilson_flow.c", "'Z' in Zeta must be defined on the whole lattice");
  error(U->type != &glattice, 1, "wilson_flow.c", "'U' in Zeta must be defined on the whole lattice");
  suNg staple, tmp1, tmp2;

  _MASTER_FOR(&glattice, i)
  {
    for (int mu = 0; mu < 4; ++mu)
    {
      _suNg_zero(staple);
      for (int nu = (mu + 1) % 4; nu != mu; nu = (nu + 1) % 4)
      {
        suNg *u1 = _4FIELD_AT(U, iup(i, mu), nu);
        suNg *u2 = _4FIELD_AT(U, iup(i, nu), mu);
        suNg *u3 = _4FIELD_AT(U, i, nu);
        _suNg_times_suNg_dagger(tmp1, *u1, *u2);
        _suNg_times_suNg_dagger(tmp2, tmp1, *u3);
#ifdef PLAQ_WEIGHTS
        _suNg_mul(tmp2, wf_plaq_weight[i * 16 + nu * 4 + mu], tmp2);
#endif
        _suNg_add_assign(staple, tmp2);

        int j = idn(i, nu);
        u1 = _4FIELD_AT(U, iup(j, mu), nu);
        u2 = _4FIELD_AT(U, j, mu);
        u3 = _4FIELD_AT(U, j, nu);
        _suNg_times_suNg(tmp1, *u2, *u1);
        _suNg_dagger_times_suNg(tmp2, tmp1, *u3);

#ifdef PLAQ_WEIGHTS
        _suNg_mul(tmp2, wf_plaq_weight[j * 16 + nu * 4 + mu], tmp2);
#endif
        _suNg_add_assign(staple, tmp2);
      }

      _suNg_times_suNg(tmp1, *_4FIELD_AT(U, i, mu), staple);

      _suNg_dagger(tmp2, tmp1);
      _suNg_sub_assign(tmp1, tmp2);
#if !defined(GAUGE_SON) && !defined(WITH_QUATERNIONS)
      double imtr;
      _suNg_trace_im(imtr, tmp1);
      imtr = imtr / NG;
      for (int k = 0; k < NG * NG; k += NG + 1)
      {
        tmp1.c[k] -= I * imtr;
      }
#endif
      _suNg_mul(tmp1, -alpha / 2., tmp1);
      _suNg_add_assign(*_4FIELD_AT(Z, i, mu), tmp1);
    }

  } //Master_for
}

void WilsonFlow1(suNg_field *V, const double epsilon)
{

  _MASTER_FOR(&glattice, ix)
  {
    for (int mu = 0; mu < 4; ++mu)
    {
      _suNg_zero(*_4FIELD_AT(ws_gf, ix, mu));
    }
  }

  Zeta(ws_gf, V, epsilon);

  _MASTER_FOR(&glattice, ix)
  {
    suNg utmp[2];
    for (int mu = 0; mu < 4; ++mu)
    {
      WF_Exp(&utmp[0], _4FIELD_AT(ws_gf, ix, mu));
      _suNg_times_suNg(utmp[1], utmp[0], *_4FIELD_AT(V, ix, mu));
      *_4FIELD_AT(V, ix, mu) = utmp[1];
    }
  }

  start_gf_sendrecv(V);
  complete_gf_sendrecv(V);
#ifdef ROTATED_SF
  apply_BCs_on_fundamental_gauge_field();
#endif
}

//define distance between complex matrices
double max_distance(suNg_field *V, suNg_field *Vprime)
{

  double d, tmp;
  suNg diff;
  _suNg_zero(diff);

  tmp = 0;
  _MASTER_FOR(&glattice, ix)
  {
    d = 0.;

    for (int mu = 0; mu < 4; ++mu)
    {
      _suNg_mul(diff, 1., *_4FIELD_AT(V, ix, mu));
      _suNg_sub_assign(diff, *_4FIELD_AT(Vprime, ix, mu));
      _suNg_sqnorm(d, diff);
      if (d > tmp)
        tmp = d;
    }
  }
  global_max(&tmp, 1);

  return tmp / (double)NG;
}

// following 1301.4388
double WilsonFlow3_adaptative(suNg_field *V, double epsilon, double delta)
{

  suNg_field_copy(u_gauge_backup, u_gauge);
  double varepsilon, d;

  _MASTER_FOR(&glattice, ix)
  {
    for (int mu = 0; mu < 4; ++mu)
    {
      _suNg_zero(*_4FIELD_AT(ws_gf, ix, mu));
      _suNg_zero(*_4FIELD_AT(ws_gf_tmp, ix, mu));
    }
  }

  start_gf_sendrecv(V);
  complete_gf_sendrecv(V);
#ifdef ROTATED_SF
  apply_BCs_on_fundamental_gauge_field();
#endif

  Zeta(ws_gf, V, epsilon / 4.); //ws_gf = Z0/4

  _MASTER_FOR(&glattice, ix)
  {
    suNg utmp[2];
    for (int mu = 0; mu < 4; ++mu)
    {
      WF_Exp(&utmp[0], _4FIELD_AT(ws_gf, ix, mu));
      _suNg_times_suNg(utmp[1], utmp[0], *_4FIELD_AT(V, ix, mu));
      *_4FIELD_AT(V, ix, mu) = utmp[1];                                             // V = exp(Z0/4) W0
      _suNg_mul(*_4FIELD_AT(ws_gf_tmp, ix, mu), -4., *_4FIELD_AT(ws_gf, ix, mu));   //ws_gf_tmp = -Z0
      _suNg_mul(*_4FIELD_AT(ws_gf, ix, mu), -17. / 9., *_4FIELD_AT(ws_gf, ix, mu)); //ws_gf =  -17*Z0/36
    }
  }

  start_gf_sendrecv(V);
  complete_gf_sendrecv(V);
#ifdef ROTATED_SF
  apply_BCs_on_fundamental_gauge_field();
#endif

  Zeta(ws_gf, V, 8. * epsilon / 9.); // ws_gf = 8 Z1 /9 - 17 Z0/36
  Zeta(ws_gf_tmp, V, 2. * epsilon);  // ws_gf_tmp = 2 Z1 - Z0

  _MASTER_FOR(&glattice, ix)
  {
    suNg utmp[4];
    for (int mu = 0; mu < 4; ++mu)
    {
      WF_Exp(&utmp[0], _4FIELD_AT(ws_gf, ix, mu));
      WF_Exp(&utmp[2], _4FIELD_AT(ws_gf_tmp, ix, mu));
      _suNg_times_suNg(utmp[1], utmp[0], *_4FIELD_AT(V, ix, mu)); // utmp[1] = exp(8 Z1/9 - 17 Z0/36) W1
      _suNg_times_suNg(utmp[3], utmp[2], *_4FIELD_AT(V, ix, mu)); // utmp[4] = exp( Z1 -  Z0) W1
      *_4FIELD_AT(V, ix, mu) = utmp[1];
      *_4FIELD_AT(Vprime, ix, mu) = utmp[3];
      _suNg_mul(*_4FIELD_AT(ws_gf, ix, mu), -1., *_4FIELD_AT(ws_gf, ix, mu));
    }
  }

  start_gf_sendrecv(V);
  complete_gf_sendrecv(V);

  start_gf_sendrecv(Vprime);
  complete_gf_sendrecv(Vprime);

#ifdef ROTATED_SF
  apply_BCs_on_fundamental_gauge_field();
#endif

  Zeta(ws_gf, V, 3. * epsilon / 4.);

  _MASTER_FOR(&glattice, ix)
  {
    suNg utmp[2];
    for (int mu = 0; mu < 4; ++mu)
    {
      WF_Exp(&utmp[0], _4FIELD_AT(ws_gf, ix, mu));
      _suNg_times_suNg(utmp[1], utmp[0], *_4FIELD_AT(V, ix, mu));
      *_4FIELD_AT(V, ix, mu) = utmp[1];
    }
  }

  start_gf_sendrecv(V);
  complete_gf_sendrecv(V);
#ifdef ROTATED_SF
  apply_BCs_on_fundamental_gauge_field();
#endif

  // now need to get the maximum of the distance
  d = max_distance(V, Vprime);

  // compute new value of epsilon
  varepsilon = 0.95 * epsilon * pow(delta / d, 1. / 3.);

  if (varepsilon < 0.1)
    epsilon = varepsilon;
  else
    epsilon = 0.1;

  if (d > delta)
  {
    suNg_field_copy(u_gauge, u_gauge_backup);
    epsilon = -1.;
    lprintf("WARNING", 0, "d > delta ! Epsilon is set to -1 in order to repeat the calculation \n");
  }

  return epsilon;
}

void WilsonFlow3(suNg_field *V, const double epsilon)
{

  _MASTER_FOR(&glattice, ix)
  {
    for (int mu = 0; mu < 4; ++mu)
    {
      _suNg_zero(*_4FIELD_AT(ws_gf, ix, mu));
    }
  }

  start_gf_sendrecv(V);
  complete_gf_sendrecv(V);
#ifdef ROTATED_SF
  apply_BCs_on_fundamental_gauge_field();
#endif

  Zeta(ws_gf, V, epsilon / 4.);

  _MASTER_FOR(&glattice, ix)
  {
    suNg utmp[2];
    for (int mu = 0; mu < 4; ++mu)
    {
      WF_Exp(&utmp[0], _4FIELD_AT(ws_gf, ix, mu));
      _suNg_times_suNg(utmp[1], utmp[0], *_4FIELD_AT(V, ix, mu));
      *_4FIELD_AT(V, ix, mu) = utmp[1];
      _suNg_mul(*_4FIELD_AT(ws_gf, ix, mu), -17. / 9., *_4FIELD_AT(ws_gf, ix, mu));
    }
  }

  start_gf_sendrecv(V);
  complete_gf_sendrecv(V);
#ifdef ROTATED_SF
  apply_BCs_on_fundamental_gauge_field();
#endif

  Zeta(ws_gf, V, 8. * epsilon / 9.);

  _MASTER_FOR(&glattice, ix)
  {
    suNg utmp[2];
    for (int mu = 0; mu < 4; ++mu)
    {
      WF_Exp(&utmp[0], _4FIELD_AT(ws_gf, ix, mu));
      _suNg_times_suNg(utmp[1], utmp[0], *_4FIELD_AT(V, ix, mu));
      *_4FIELD_AT(V, ix, mu) = utmp[1];
      _suNg_mul(*_4FIELD_AT(ws_gf, ix, mu), -1., *_4FIELD_AT(ws_gf, ix, mu));
    }
  }

  start_gf_sendrecv(V);
  complete_gf_sendrecv(V);
#ifdef ROTATED_SF
  apply_BCs_on_fundamental_gauge_field();
#endif

  Zeta(ws_gf, V, 3. * epsilon / 4.);

  _MASTER_FOR(&glattice, ix)
  {
    suNg utmp[2];
    for (int mu = 0; mu < 4; ++mu)
    {
      WF_Exp(&utmp[0], _4FIELD_AT(ws_gf, ix, mu));
      _suNg_times_suNg(utmp[1], utmp[0], *_4FIELD_AT(V, ix, mu));
      *_4FIELD_AT(V, ix, mu) = utmp[1];
    }
  }

  start_gf_sendrecv(V);
  complete_gf_sendrecv(V);
#ifdef ROTATED_SF
  apply_BCs_on_fundamental_gauge_field();
#endif
}

static void WF_plaq(double *ret, suNg_field *V, int ix, int mu, int nu)
{
  int iy, iz;
  suNg *v1, *v2, *v3, *v4, w1, w2, w3;

  iy = iup(ix, mu);
  iz = iup(ix, nu);

  v1 = _4FIELD_AT(V, ix, mu);
  v2 = _4FIELD_AT(V, iy, nu);
  v3 = _4FIELD_AT(V, iz, mu);
  v4 = _4FIELD_AT(V, ix, nu);

  _suNg_times_suNg(w1, (*v1), (*v2));
  _suNg_times_suNg(w2, (*v4), (*v3));
  _suNg_times_suNg_dagger(w3, w1, w2);

  _suNg_trace_re(*ret, w3);

#ifdef PLAQ_WEIGHTS
  *ret *= wf_plaq_weight[ix * 16 + nu * 4 + mu];
#endif
}

double WF_E(suNg_field *V)
{
  double E = 0.;

  _MASTER_FOR_SUM(&glattice, ix, E)
  {
    double p;
    for (int mu = 0; mu < 4; mu++)
      for (int nu = mu + 1; nu < 4; nu++)
      {
        WF_plaq(&p, V, ix, mu, nu);
        E += ((double)NG) - p;
      }
  }

  E *= 2. / ((double)GLB_VOLUME);

  global_sum(&E, 1);

  return E;
}

void WF_E_T(double *E, suNg_field *V)
{
  int gt, t, x, y, z, ix;
  int mu, nu;
  double p;

  for (t = 0; t < 2 * GLB_T; t++)
    E[t] = 0.;

  for (t = 0; t < T; t++)
  {
    gt = t + zerocoord[0];
    for (x = 0; x < X; x++)
      for (y = 0; y < Y; y++)
        for (z = 0; z < Z; z++)
        {
          mu = 0;
          ix = ipt(t, x, y, z);
          for (nu = 1; nu < 4; nu++)
          {
            WF_plaq(&p, V, ix, mu, nu);
#ifdef PLAQ_WEIGHTS
            E[2 * gt] += ((double)(NG)) * wf_plaq_weight[ix * 16 + nu * 4 + mu] - p;
#else
            E[2 * gt] += NG - p;
#endif
          }
          for (mu = 1; mu < 3; mu++)
            for (nu = mu + 1; nu < 4; nu++)
            {
              WF_plaq(&p, V, ix, mu, nu);
#ifdef PLAQ_WEIGHTS
              E[2 * gt + 1] += ((double)(NG)) * wf_plaq_weight[ix * 16 + nu * 4 + mu] - p;
#else
              E[2 * gt + 1] += NG - p;
#endif
            }
        }
    E[2 * gt] /= (3. * GLB_VOL3);
    E[2 * gt + 1] /= (3. * GLB_VOL3);
  }

  global_sum(E, 2 * GLB_T);
}

/* This gives F_{\mu\nu}^A */
static void WF_clover_F(suNg_algebra_vector *F, suNg_field *V, int ix, int mu, int nu)
{
  int iy, iz, iw;
  suNg *v1, *v2, *v3, *v4, w1, w2, w3;

  _suNg_unit(w3);
  _suNg_mul(w3, -4., w3);

  iy = iup(ix, mu);
  iz = iup(ix, nu);

  v1 = _4FIELD_AT(V, ix, mu);
  v2 = _4FIELD_AT(V, iy, nu);
  v3 = _4FIELD_AT(V, iz, mu);
  v4 = _4FIELD_AT(V, ix, nu);

  _suNg_times_suNg(w1, (*v1), (*v2));
  _suNg_times_suNg_dagger(w2, w1, (*v3));
  _suNg_times_suNg_dagger(w1, w2, (*v4));
#ifdef PLAQ_WEIGHTS
  _suNg_mul(w1, wf_plaq_weight[ix * 16 + nu * 4 + mu], w1);
#endif
  _suNg_add_assign(w3, w1);

  iy = idn(ix, mu);
  iz = iup(iy, nu);

  v1 = _4FIELD_AT(V, ix, nu);
  v2 = _4FIELD_AT(V, iz, mu);
  v3 = _4FIELD_AT(V, iy, nu);
  v4 = _4FIELD_AT(V, iy, mu);

  _suNg_times_suNg_dagger(w1, (*v1), (*v2));
  _suNg_times_suNg_dagger(w2, w1, (*v3));
  _suNg_times_suNg(w1, w2, (*v4));
#ifdef PLAQ_WEIGHTS
  _suNg_mul(w1, wf_plaq_weight[iy * 16 + nu * 4 + mu], w1);
#endif
  _suNg_add_assign(w3, w1);

  iy = idn(ix, mu);
  iz = idn(iy, nu);
  iw = idn(ix, nu);

  v1 = _4FIELD_AT(V, iy, mu);
  v2 = _4FIELD_AT(V, iz, nu);
  v3 = _4FIELD_AT(V, iz, mu);
  v4 = _4FIELD_AT(V, iw, nu);

  _suNg_times_suNg(w1, (*v2), (*v1));
  _suNg_dagger_times_suNg(w2, w1, (*v3));
  _suNg_times_suNg(w1, w2, (*v4));
#ifdef PLAQ_WEIGHTS
  _suNg_mul(w1, wf_plaq_weight[iz * 16 + nu * 4 + mu], w1);

#endif
  _suNg_add_assign(w3, w1);

  iy = idn(ix, nu);
  iz = iup(iy, mu);

  v1 = _4FIELD_AT(V, iy, nu);
  v2 = _4FIELD_AT(V, iy, mu);
  v3 = _4FIELD_AT(V, iz, nu);
  v4 = _4FIELD_AT(V, ix, mu);

  _suNg_dagger_times_suNg(w1, (*v1), (*v2));
  _suNg_times_suNg(w2, w1, (*v3));
  _suNg_times_suNg_dagger(w1, w2, (*v4));
#ifdef PLAQ_WEIGHTS
  _suNg_mul(w1, wf_plaq_weight[iy * 16 + nu * 4 + mu], w1);
#endif
  _suNg_add_assign(w3, w1);

  _fund_algebra_project(*F, w3);

  _algebra_vector_mul_g(*F, 1 / 4., *F);
}

double WF_Esym(suNg_field *V)
{
  double E = 0.;

  _MASTER_FOR_SUM(&glattice, ix, E)
  {
    suNg_algebra_vector clover;
    double p;
    for (int mu = 0; mu < 4; mu++)
      for (int nu = mu + 1; nu < 4; nu++)
      {
        WF_clover_F(&clover, V, ix, mu, nu);
        _algebra_vector_sqnorm_g(p, clover);
        E += p;
      }
  }

  E *= _FUND_NORM2 / ((double)GLB_VOLUME);

  global_sum(&E, 1);
  return E;
}

void WF_Esym_T(double *E, suNg_field *V)
{
  int gt, t, x, y, z, ix;
  int mu, nu;
  suNg_algebra_vector clover;
  double p;

  for (t = 0; t < 2 * GLB_T; t++)
    E[t] = 0.;

  for (t = 0; t < T; t++)
  {
    gt = t + zerocoord[0];
    for (x = 0; x < X; x++)
      for (y = 0; y < Y; y++)
        for (z = 0; z < Z; z++)
        {
          mu = 0;
          ix = ipt(t, x, y, z);
          for (nu = 1; nu < 4; nu++)
          {
            WF_clover_F(&clover, V, ix, mu, nu);
            _algebra_vector_sqnorm_g(p, clover);
            E[2 * gt] += p;
          }

          for (mu = 1; mu < 4; mu++)
            for (nu = mu + 1; nu < 4; nu++)
            {
              WF_clover_F(&clover, V, ix, mu, nu);
              _algebra_vector_sqnorm_g(p, clover);
              E[2 * gt + 1] += p;
            }
        }
    E[2 * gt] *= _FUND_NORM2 / (6. * GLB_VOL3);
    E[2 * gt + 1] *= _FUND_NORM2 / (6. * GLB_VOL3);
  }

  global_sum(E, 2 * GLB_T);
}

/*
  q = 1/(16 \pi^2) \epsilon_{\mu\nu\rho\sigma} \tr F_{\mu\nu} F_{\rho\sigma}
*/
double WF_topo(suNg_field *V)
{
  double TC = 0.;

  _MASTER_FOR_SUM(&glattice, ix, TC)
  {
    suNg_algebra_vector F1, F2;
    WF_clover_F(&F1, V, ix, 1, 2);
    WF_clover_F(&F2, V, ix, 0, 3);
    for (int i = 0; i < NG * NG - 1; i++)
      TC += F1.c[i] * F2.c[i];

    WF_clover_F(&F1, V, ix, 1, 3);
    WF_clover_F(&F2, V, ix, 0, 2);
    for (int i = 0; i < NG * NG - 1; i++)
      TC -= F1.c[i] * F2.c[i];

    WF_clover_F(&F1, V, ix, 0, 1);
    WF_clover_F(&F2, V, ix, 2, 3);
    for (int i = 0; i < NG * NG - 1; i++)
      TC += F1.c[i] * F2.c[i];
  }

  TC *= _FUND_NORM2 / (4. * M_PI * M_PI);

  global_sum(&TC, 1);

  return TC;
}

void WF_adaptive_full_measure(suNg_field *V, double *tmax, double *eps, double *delta, int nmeas)
{

  double TC;
#if defined(BC_T_ANTIPERIODIC) || defined(BC_T_PERIODIC) && !defined(PURE_GAUGE_ANISOTROPY)
  double E, Esym;
#else
  int j;
  double E[2 * GLB_T];
  double Esym[2 * GLB_T];
  double Eavg[2];
  double Esymavg[2];
#endif

  int k;
  double t = 0.;
  double epsilon = *eps;
  double dt = *tmax / nmeas;

  TC = WF_topo(V);

#if defined(BC_T_ANTIPERIODIC) || defined(BC_T_PERIODIC) && !defined(PURE_GAUGE_ANISOTROPY)
  E = WF_E(V);
  Esym = WF_Esym(V);
  lprintf("WILSONFLOW", 0, "WF (t,E,t2*E,Esym,t2*Esym,TC) = %e %e %e %e %e %e\n", t, E, t * t * E, Esym, t * t * Esym, TC);
#else

  WF_E_T(E, V);
  WF_Esym_T(Esym, V);
  Eavg[0] = Eavg[1] = Esymavg[0] = Esymavg[1] = 0.0;
  for (j = 0; j < GLB_T; j++)
  {
    lprintf("WILSONFLOW", 0, "WF (T,t,Etime,Espace,Esymtime,Esymspace) = %d %e %e %e %e %e\n", j, t, E[2 * j], E[2 * j + 1], Esym[2 * j], Esym[2 * j + 1]);
    Eavg[0] += E[2 * j];
    Eavg[1] += E[2 * j + 1];
    Esymavg[0] += Esym[2 * j];
    Esymavg[1] += Esym[2 * j + 1];
  }
  Eavg[0] /= GLB_T - 2;
  Eavg[1] /= GLB_T - 3;
  Esymavg[0] /= GLB_T - 2;
  Esymavg[1] /= GLB_T - 3;

  lprintf("WILSONFLOW", 0, "WF avg (t,Etime,Espace,Esymtime,Esymspace,Pltime,Plspace,TC) = %e %e %e %e %e %e %e %e\n", t, Eavg[0], Eavg[1], Esymavg[0], Esymavg[1], (NG - Eavg[0]), (NG - Eavg[1]), TC);

#endif
  k = 1;
  double epsilon_new = 0;
  while (t < *tmax)
  {
    if (t + epsilon > (double)k * dt)
      epsilon = (double)k * dt - t;

    epsilon_new = WilsonFlow3_adaptative(V, epsilon, *delta);

    if (fabs(epsilon_new + 1.) > 1e-7)
      t = t + epsilon;

    if (fabs(t - (double)k * dt) < 1e-7)
    {
      k = k + 1;

      TC = WF_topo(V);

#if defined(BC_T_ANTIPERIODIC) || defined(BC_T_PERIODIC) && !defined(PURE_GAUGE_ANISOTROPY)

      E = WF_E(V);
      Esym = WF_Esym(V);
      lprintf("WILSONFLOW", 0, "WF (t,E,t2*E,Esym,t2*Esym,TC) = %e %e %e %e %e %e\n", t, E, t * t * E, Esym, t * t * Esym, TC);
#else

      WF_E_T(E, V);
      WF_Esym_T(Esym, V);
      Eavg[0] = Eavg[1] = Esymavg[0] = Esymavg[1] = 0.0;
      for (j = 1; j < GLB_T - 1; j++)
      {
        lprintf("WILSONFLOW", 0, "WF (T,t,Etime,Espace,Esymtime,Esymspace) = %d %e %e %e %e %e\n", j, t, E[2 * j], E[2 * j + 1], Esym[2 * j], Esym[2 * j + 1]);
        Eavg[0] += E[2 * j];
        Eavg[1] += E[2 * j + 1];
        Esymavg[0] += Esym[2 * j];
        Esymavg[1] += Esym[2 * j + 1];
      }

      Eavg[0] /= GLB_T - 2;
      Eavg[1] /= GLB_T - 3;
      Esymavg[0] /= GLB_T - 2;
      Esymavg[1] /= GLB_T - 3;

      lprintf("WILSONFLOW", 0, "WF avg (t,Etime,Espace,Esymtime,Esymspace,Pltime,Plspace,TC) = %e %e %e %e %e %e %e %e\n", t, Eavg[0], Eavg[1], Esymavg[0], Esymavg[1], (NG - Eavg[0]), (NG - Eavg[1]), TC);

#endif
    }
    if (fabs(epsilon_new + 1.) > 1e-7)
      epsilon = epsilon_new;
    if (fabs(epsilon_new + 1.) < 1e-7)
      epsilon = epsilon / 2;
  }
}