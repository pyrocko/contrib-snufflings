/*
 * Copyright (c) 2010 Gertjan van Zwieten
 * Copyright (c) 2010 Delft University of Technology, The Netherlands
 *
 * This file is part of DeForM, the Delft Forward Modeling toolkit.
 *
 * DeForM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * DeForM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Publications that contain results produced by the DeForM software should
 * contain an acknowledgment. For example: The forward modeling was performed
 * using the freely available DeForM software developed by the Delft Institute
 * of Earth Observation and Space Systems (DEOS), Delft University of
 * Technology, or include a reference to Gertjan van Zwieten.
 *
 */

#ifdef DLL
  #define EXPORT __declspec(dllexport)
#else
  #define EXPORT
#endif

#define _USE_MATH_DEFINES

#include <math.h>
#include <stdio.h>

__inline int isTiny( double x ) {
  return x < 1e-10 && x > -1e-10;
}

// +-----------------------+
// |                       |
// |      ANGLE TYPE       |
// |                       |
// +-----------------------+

typedef struct { // Angle
  double sin, cos;
} Angle;
__inline Angle from_deg( const double deg ) {
  const double rad = deg * M_PI / 180.;
  Angle angle = { sin( rad ), cos( rad ) };
  return angle;
}

// +-----------------------+
// |                       |
// |      VECTOR TYPE      |
// |                       |
// +-----------------------+

typedef struct { // Vector
  double x, y, z;
} Vector;
__inline void Vector_copy( Vector *self, const Vector other ) {
  self->x = other.x;
  self->y = other.y;
  self->z = other.z;
}
__inline void Vector_flip( Vector *self ) {
  self->x = -self->x;
  self->y = -self->y;
  self->z = -self->z;
}
__inline void Vector_iadd( Vector *self, const Vector other ) {
  self->x += other.x;
  self->y += other.y;
  self->z += other.z;
}
__inline void Vector_isub( Vector *self, const Vector other ) {
  self->x -= other.x;
  self->y -= other.y;
  self->z -= other.z;
}
__inline void Vector_imul( Vector *self, const double factor ) {
  self->x *= factor;
  self->y *= factor;
  self->z *= factor;
}
__inline void Vector_irot_yz( Vector *self, const Angle angle ) {
  const double y = self->y, z = self->z;
  self->y = y * angle.cos - z * angle.sin;
  self->z = z * angle.cos + y * angle.sin;
}
__inline void Vector_irot_zy( Vector *self, const Angle angle ) {
  const double y = self->y, z = self->z;
  self->y = y * angle.cos + z * angle.sin;
  self->z = z * angle.cos - y * angle.sin;
}
__inline void Vector_irot_zx( Vector *self, const Angle angle ) {
  const double z = self->z, x = self->x;
  self->z = z * angle.cos - x * angle.sin;
  self->x = x * angle.cos + z * angle.sin;
}
__inline void Vector_irot_xz( Vector *self, const Angle angle ) {
  const double z = self->z, x = self->x;
  self->z = z * angle.cos + x * angle.sin;
  self->x = x * angle.cos - z * angle.sin;
}
__inline void Vector_irot_xy( Vector *self, const Angle angle ) {
  const double x = self->x, y = self->y;
  self->x = x * angle.cos - y * angle.sin;
  self->y = y * angle.cos + x * angle.sin;
}
__inline void Vector_irot_yx( Vector *self, const Angle angle ) {
  const double x = self->x, y = self->y;
  self->x = x * angle.cos + y * angle.sin;
  self->y = y * angle.cos - x * angle.sin;
}
__inline Vector Vector_mul( Vector self, double factor ) {
  Vector retval = { self.x * factor, self.y * factor, self.z * factor };
  return retval;
}

// +-----------------------+
// |                       |
// |      MATRIX TYPE      |
// |                       |
// +-----------------------+

typedef struct { // Matrix
  Vector x, y, z;
} Matrix;
__inline void Matrix_iadd( Matrix *self, const Matrix other ) {
  Vector_iadd( &(self->x), other.x );
  Vector_iadd( &(self->y), other.y );
  Vector_iadd( &(self->z), other.z );
}
__inline void Matrix_isub( Matrix *self, const Matrix other ) {
  Vector_isub( &(self->x), other.x );
  Vector_isub( &(self->y), other.y );
  Vector_isub( &(self->z), other.z );
}
__inline void Matrix_imul( Matrix *self, const double factor ) {
  Vector_imul( &(self->x), factor );
  Vector_imul( &(self->y), factor );
  Vector_imul( &(self->z), factor );
}
__inline void Matrix_copy( Matrix *self, const Matrix other ) {
  Vector_copy( &(self->x), other.x );
  Vector_copy( &(self->y), other.y );
  Vector_copy( &(self->z), other.z );
}
__inline void Matrix_irot_xy( Matrix *self, const Angle angle ) {
  const double yyxx = ( self->y.y - self->x.x ) * angle.sin;
  const double xyyx = ( self->x.y + self->y.x ) * angle.sin;
  const double diag =  yyxx * angle.sin + xyyx * angle.cos;
  const double offd = -xyyx * angle.sin + yyxx * angle.cos;
  const double xz = self->x.z;
  const double yz = self->y.z;
  self->x.x += diag;
  self->x.y += offd;
  self->x.z = xz * angle.cos + yz * angle.sin;
  self->y.x += offd;
  self->y.y -= diag;
  self->y.z = yz * angle.cos - xz * angle.sin;
  Vector_irot_yx( &(self->z), angle );
}

// +-----------------------+
// |                       |
// |   HELPER FUNCTIONS    |
// |                       |
// +-----------------------+

typedef struct { // OkadaConsts
  double sin, cos;
  double z, d, p, q, h;
  double xi, eta;
  double ytilde, dtilde, ctilde;
  double R, R2, R3, R5, X, theta;
  double X11, Y11, X32, Y32, Z32, X53, Y53, Z53;
  double Y0, Z0;
  double I1, I2, I3, I4;
  // more
  double D11, K1, K2, K3, K4;
  double J1, J2, J3, J4, J5, J6;
  double E, F, G, H, P, Q;
  double Eprime, Fprime, Gprime, Hprime, Pprime, Qprime;
  // elasticity
  double ALPHA;
} OkadaConsts;
__inline void okada_get_consts( OkadaConsts *consts, const Vector coord, const double c, const Angle dip, const double length, const double width, const int n, const double poisson ) {
  consts->sin    = dip.sin;
  consts->cos    = dip.cos;
  consts->z      = coord.z;
  consts->d      = c - consts->z;
  consts->p      = coord.y * dip.cos + consts->d * dip.sin;
  consts->q      = coord.y * dip.sin - consts->d * dip.cos;
  consts->h      = consts->q * dip.cos - consts->z;
  consts->xi     = coord.x   - ( n & 1 ? 0.5 : -0.5 ) * length;
  consts->eta    = consts->p - ( n & 2 ? 1.0 :  0.0 ) * width;
  consts->ytilde = consts->eta * dip.cos + consts->q * dip.sin;
  consts->dtilde = consts->eta * dip.sin - consts->q * dip.cos;
  consts->ctilde = consts->dtilde + consts->z;
  consts->R2     = consts->xi * consts->xi + consts->eta * consts->eta + consts->q * consts->q;
  consts->R      = sqrt( consts->R2 );
  consts->R3     = consts->R * consts->R2;
  consts->R5     = consts->R2 * consts->R3;
  consts->X      = hypot( consts->xi, consts->q );
  consts->theta  = atan( ( consts->xi * consts->eta ) / ( consts->q * consts->R ) );
  consts->X11    = 1. / ( consts->R * ( consts->R + consts->xi ) );
  consts->Y11    = 1. / ( consts->R * ( consts->R + consts->eta ) );
  consts->X32    = ( 2. * consts->R + consts->xi  ) / ( consts->R3 * ( consts->R + consts->xi  ) * ( consts->R + consts->xi ) );
  consts->Y32    = ( 2. * consts->R + consts->eta ) / ( consts->R3 * ( consts->R + consts->eta ) * ( consts->R + consts->eta ) );
  consts->Z32    = dip.sin / consts->R3 - consts->h * consts->Y32;
  consts->X53    = ( 8. * consts->R2 + 9. * consts->R * consts->xi  + 3. * consts->xi  * consts->xi  ) / ( consts->R5 * ( consts->R + consts->xi  ) * ( consts->R + consts->xi  ) * ( consts->R + consts->xi  ) );
  consts->Y53    = ( 8. * consts->R2 + 9. * consts->R * consts->eta + 3. * consts->eta * consts->eta ) / ( consts->R5 * ( consts->R + consts->eta ) * ( consts->R + consts->eta ) * ( consts->R + consts->eta ) );
  consts->Z53    = 3. * dip.sin / consts->R5 - consts->h * consts->Y53;
  consts->Y0     = consts->Y11 - consts->xi * consts->xi * consts->Y32;
  consts->Z0     = consts->Z32 - consts->xi * consts->xi * consts->Z53;
  if ( isTiny( dip.cos ) ) {
    consts->I4   = .5 * consts->xi * consts->ytilde / ( ( consts->R + consts->dtilde ) * ( consts->R + consts->dtilde ) );
    consts->I3   = .5 * consts->eta / ( consts->R + consts->dtilde ) + .5 * consts->ytilde * consts->q / ( ( consts->R + consts->dtilde ) * ( consts->R + consts->dtilde ) ) - .5 * log( consts->R + consts->eta );
  }
  else {
    consts->I4   = consts->xi * dip.sin / dip.cos / ( consts->R + consts->dtilde ) + 2. * atan( ( consts->eta * ( consts->X + consts->q * dip.cos ) + consts->X * ( consts->R + consts->X ) * dip.sin ) / ( consts->xi * ( consts->R + consts->X ) * dip.cos ) ) / ( dip.cos * dip.cos );
    consts->I3   = consts->ytilde / dip.cos / ( consts->R + consts->dtilde ) - ( log( consts->R + consts->eta ) - dip.sin * log( consts->R + consts->dtilde ) ) / ( dip.cos * dip.cos );
  }
  consts->I2     = log( consts->R + consts->dtilde ) + consts->I3 * dip.sin;
  consts->I1     = -consts->xi * dip.cos / ( consts->R + consts->dtilde ) - consts->I4 * dip.sin;
  consts->ALPHA  = .5 / ( 1. - poisson );
};
__inline void okada_get_more_consts( OkadaConsts *consts, const Vector coord, const double c, const Angle dip, const double length, const double width, const int n, const double poisson ) {
  okada_get_consts( consts, coord, c, dip, length, width, n, poisson );
  consts->D11    = 1. / ( consts->R * ( consts->R + consts->dtilde ) );
  if ( isTiny( dip.cos ) ) {
    consts->K1   = consts->xi * consts->q / ( consts->R + consts->dtilde ) * consts->D11;
    consts->K3   = dip.sin / ( consts->R + consts->dtilde ) * ( consts->xi * consts->xi * consts->D11 - 1. );
  }
  else {
    consts->K1   = consts->xi * ( consts->D11 - consts->Y11 * dip.sin ) / dip.cos;
    consts->K3   = ( consts->q * consts->Y11  - consts->ytilde * consts->D11 ) / dip.cos;
  }
  consts->K2     = 1. / consts->R + consts->K3 * dip.sin;
  consts->K4     = consts->xi * consts->Y11 * dip.cos - consts->K1 * dip.sin;
  consts->J2     = consts->xi * consts->ytilde / ( consts->R + consts->dtilde ) * consts->D11;
  consts->J5     = -( consts->dtilde + consts->ytilde * consts->ytilde / ( consts->R + consts->dtilde ) ) * consts->D11;
  if ( isTiny( dip.cos ) ) {
    consts->J3   = -consts->xi * ( consts->q * consts->q * consts->D11 - .5 ) / ( ( consts->R + consts->dtilde ) * ( consts->R + consts->dtilde ) );
    consts->J6   = -consts->ytilde * ( consts->xi * consts->xi * consts->D11 - .5 ) / ( ( consts->R + consts->dtilde ) * ( consts->R + consts->dtilde ) );
  }
  else {
    consts->J3   = ( consts->K1 - consts->J2 * dip.sin ) / dip.cos;
    consts->J6   = ( consts->K3 - consts->J5 * dip.sin ) / dip.cos;
  }
  consts->J4     = -consts->xi * consts->Y11 - consts->J2 * dip.cos + consts->J3 * dip.sin;
  consts->J1     = consts->J5 * dip.cos - consts->J6 * dip.sin;
  consts->E      = dip.sin / consts->R - consts->ytilde * consts->q / consts->R3;
  consts->F      = consts->dtilde / consts->R3 + consts->xi * consts->xi * consts->Y32 * dip.sin;
  consts->G      = 2. * consts->X11 * dip.sin - consts->ytilde * consts->q * consts->X32;
  consts->H      = consts->dtilde * consts->q * consts->X32 + consts->xi * consts->q * consts->Y32 * dip.sin;
  consts->P      = dip.cos / consts->R3 + consts->q * consts->Y32 * dip.sin;
  consts->Q      = 3. * consts->ctilde * consts->dtilde / consts->R5 - ( consts->z * consts->Y32 + consts->Z32 + consts->Z0 ) * dip.sin;
  consts->Eprime = dip.cos / consts->R + consts->dtilde * consts->q / consts->R3;
  consts->Fprime = consts->ytilde / consts->R3 + consts->xi * consts->xi * consts->Y32 * dip.cos;
  consts->Gprime = 2. * consts->X11 * dip.cos + consts->dtilde * consts->q * consts->X32;
  consts->Hprime = consts->ytilde * consts->q * consts->X32 + consts->xi * consts->q * consts->Y32 * dip.cos;
  consts->Pprime = dip.sin / consts->R3 - consts->q * consts->Y32 * dip.cos;
  consts->Qprime = 3. * consts->ctilde * consts->ytilde / consts->R5 + consts->q * consts->Y32 - ( consts->z * consts->Y32 + consts->Z32 + consts->Z0 ) * dip.cos;
};
__inline void okada_add_disp_A( Vector *disp, const OkadaConsts *consts, const Vector U ) {
  if ( U.x ) {
    disp->x += U.x * ( .5 * consts->theta + .5 * consts->ALPHA * consts->xi * consts->q * consts->Y11 );
    disp->y += U.x * ( .5 * consts->ALPHA * consts->q / consts->R );
    disp->z += U.x * ( .5 * ( 1.0 - consts->ALPHA ) * log( consts->R + consts->eta ) - .5 * consts->ALPHA * consts->q * consts->q * consts->Y11 );
  }
  if ( U.y ) {
    disp->x += U.y * ( .5 * consts->ALPHA * consts->q / consts->R );
    disp->y += U.y * ( .5 * consts->theta + .5 * consts->ALPHA * consts->eta * consts->q * consts->X11 );
    disp->z += U.y * ( .5 * ( 1 - consts->ALPHA ) * log( consts->R + consts->xi ) - .5 * consts->ALPHA * consts->q * consts->q * consts->X11 );
  }
  if ( U.z ) {
    disp->x += U.z * ( -.5 * ( 1.0 - consts->ALPHA ) * log( consts->R + consts->eta ) - .5 * consts->ALPHA * consts->q * consts->q * consts->Y11 );
    disp->y += U.z * ( -.5 * ( 1.0 - consts->ALPHA ) * log( consts->R + consts->xi  ) - .5 * consts->ALPHA * consts->q * consts->q * consts->X11 );
    disp->z += U.z * ( .5 * consts->theta - .5 * consts->ALPHA * consts->q * ( consts->eta * consts->X11 + consts->xi * consts->Y11 ) );
  }
}
__inline void okada_add_disp_B( Vector *disp, const OkadaConsts *consts, const Vector U ) {
  if ( U.x ) {
    disp->x += U.x * ( -consts->xi * consts->q * consts->Y11 - consts->theta - ( 1.0 - consts->ALPHA ) / consts->ALPHA * consts->I1 * consts->sin );
    disp->y += U.x * ( -consts->q / consts->R + ( 1.0 - consts->ALPHA ) / consts->ALPHA * consts->ytilde / ( consts->R + consts->dtilde ) * consts->sin );
    disp->z += U.x * ( consts->q * consts->q * consts->Y11 - ( 1.0 - consts->ALPHA ) / consts->ALPHA * consts->I2 * consts->sin );
  }
  if ( U.y ) {
    disp->x += U.y * ( -consts->q / consts->R + ( 1.0 - consts->ALPHA ) / consts->ALPHA * consts->I3 * consts->sin * consts->cos );
    disp->y += U.y * ( -consts->eta * consts->q * consts->X11 - consts->theta - ( 1.0 - consts->ALPHA ) / consts->ALPHA * consts->xi / ( consts->R + consts->dtilde ) * consts->sin * consts->cos );
    disp->z += U.y * ( consts->q * consts->q * consts->X11 + ( 1.0 - consts->ALPHA ) / consts->ALPHA * consts->I4 * consts->sin * consts->cos );
  }
  if ( U.z ) {
    disp->x += U.z * ( consts->q * consts->q * consts->Y11 - ( 1.0 - consts->ALPHA ) / consts->ALPHA * consts->I3 * consts->sin * consts->sin );
    disp->y += U.z * ( consts->q * consts->q * consts->X11 + ( 1.0 - consts->ALPHA ) / consts->ALPHA * consts->xi / ( consts->R + consts->dtilde ) * consts->sin * consts->sin );
    disp->z += U.z * ( consts->q * ( consts->eta * consts->X11 + consts->xi * consts->Y11 ) - consts->theta - ( 1.0 - consts->ALPHA ) / consts->ALPHA * consts->I4 * consts->sin * consts->sin );
  }
}
__inline void okada_add_disp_C( Vector *disp, const OkadaConsts *consts, const Vector U ) {
  if ( U.x ) {
    disp->x += U.x * ( ( 1.0 - consts->ALPHA ) * consts->xi * consts->Y11 * consts->cos - consts->ALPHA * consts->xi * consts->q * consts->Z32 );
    disp->y += U.x * ( ( 1.0 - consts->ALPHA ) * ( consts->cos / consts->R + 2.0 * consts->q * consts->Y11 * consts->sin ) - consts->ALPHA * consts->ctilde * consts->q / consts->R3 );
    disp->z += U.x * ( ( 1.0 - consts->ALPHA ) * consts->q * consts->Y11 * consts->cos - consts->ALPHA * ( consts->ctilde * consts->eta / consts->R3 - consts->z * consts->Y11 + consts->xi * consts->xi * consts->Z32 ) );
  }
  if ( U.y ) {
    disp->x += U.y * ( ( 1.0 - consts->ALPHA ) * consts->cos / consts->R - consts->q * consts->Y11 * consts->sin - consts->ALPHA * consts->ctilde * consts->q / consts->R3 );
    disp->y += U.y * ( ( 1.0 - consts->ALPHA ) * consts->ytilde * consts->X11 - consts->ALPHA * consts->ctilde * consts->eta * consts->q * consts->X32 );
    disp->z += U.y * ( -consts->dtilde * consts->X11 - consts->xi * consts->Y11 * consts->sin - consts->ALPHA * consts->ctilde * ( consts->X11 - consts->q * consts->q * consts->X32 ) );
  }
  if ( U.z ) {
    disp->x += U.z * ( -( 1.0 - consts->ALPHA ) * ( consts->sin / consts->R + consts->q * consts->Y11 * consts->cos ) - consts->ALPHA * ( consts->z * consts->Y11 - consts->q * consts->q * consts->Z32 ) );
    disp->y += U.z * ( ( 1.0 - consts->ALPHA ) * 2.0 * consts->xi * consts->Y11 * consts->sin + consts->dtilde * consts->X11 - consts->ALPHA * consts->ctilde * ( consts->X11 - consts->q * consts->q * consts->X32 ) );
    disp->z += U.z * ( ( 1.0 - consts->ALPHA ) * ( consts->ytilde * consts->X11 + consts->xi * consts->Y11 * consts->cos ) + consts->ALPHA * consts->q * ( consts->ctilde * consts->eta * consts->X32 + consts->xi * consts->Z32 ) );
  }
}
__inline void okada_add_grad_Axy( Matrix *grad, const OkadaConsts *consts, const Vector U ) {
  if ( U.x ) {
    grad->x.x += U.x * ( -.5 * ( 1. - consts->ALPHA ) * consts->q * consts->Y11  - .5 * consts->ALPHA * consts->xi * consts->xi * consts->q * consts->Y32 );
    grad->x.y += U.x * ( -.5 * consts->ALPHA * consts->xi * consts->q / consts->R3 );
    grad->x.z += U.x * (  .5 * ( 1. - consts->ALPHA ) * consts->xi * consts->Y11 + .5 * consts->ALPHA * consts->xi * consts->q * consts->q * consts->Y32 );
    grad->y.x += U.x * (  .5 * ( 1. - consts->ALPHA ) * consts->xi * consts->Y11 * consts->sin + .5 * consts->dtilde * consts->X11 + .5 * consts->ALPHA * consts->xi * consts->F );
    grad->y.y += U.x * (  .5 * consts->ALPHA * consts->E );
    grad->y.z += U.x * (  .5 * ( 1. - consts->ALPHA ) * ( consts->cos / consts->R + consts->q * consts->Y11 * consts->sin ) - .5 * consts->ALPHA * consts->q * consts->F );
  }
  if ( U.y ) {
    grad->x.x += U.y * ( -.5 * consts->ALPHA * consts->xi * consts->q / consts->R3 );
    grad->x.y += U.y * ( -.5 * consts->q * consts->Y11 - .5 * consts->ALPHA * consts->eta * consts->q / consts->R3 );
    grad->x.z += U.y * (  .5 * ( 1. - consts->ALPHA ) / consts->R + .5 * consts->ALPHA * consts->q * consts->q / consts->R3 );
    grad->y.x += U.y * (  .5 * consts->ALPHA * consts->E );
    grad->y.y += U.y * (  .5 * ( 1. - consts->ALPHA ) * consts->dtilde * consts->X11 + .5 * consts->xi * consts->Y11 * consts->sin + .5 * consts->ALPHA * consts->eta * consts->G );
    grad->y.z += U.y * (  .5 * ( 1. - consts->ALPHA ) * consts->ytilde * consts->X11 - .5 * consts->ALPHA * consts->q * consts->G );
  }
  if ( U.z ) {
    grad->x.x += U.z * ( -.5 * ( 1. - consts->ALPHA ) * consts->xi * consts->Y11 + .5 * consts->ALPHA * consts->xi * consts->q * consts->q * consts->Y32 );
    grad->x.y += U.z * ( -.5 * ( 1. - consts->ALPHA ) / consts->R + .5 * consts->ALPHA * consts->q * consts->q / consts->R3 );
    grad->x.z += U.z * ( -.5 * ( 1. - consts->ALPHA ) * consts->q * consts->Y11 - .5 * consts->ALPHA * consts->q * consts->q * consts->q * consts->Y32 );
    grad->y.x += U.z * ( -.5 * ( 1. - consts->ALPHA ) * ( consts->cos / consts->R + consts->q * consts->Y11 * consts->sin ) - .5 * consts->ALPHA * consts->q * consts->F );
    grad->y.y += U.z * ( -.5 * ( 1. - consts->ALPHA ) * consts->ytilde * consts->X11 - .5 * consts->ALPHA * consts->q * consts->G );
    grad->y.z += U.z * (  .5 * ( 1. - consts->ALPHA ) * ( consts->dtilde * consts->X11 + consts->xi * consts->Y11 * consts->sin ) + .5 * consts->ALPHA * consts->q * consts->H );
  }
}
__inline void okada_add_grad_Az( Vector *gradz, const OkadaConsts *consts, const Vector U ) {
  if ( U.x ) {
    gradz->x += U.x * (  .5 * ( 1. - consts->ALPHA ) * consts->xi * consts->Y11 * consts->cos + .5 * consts->ytilde * consts->X11 + .5 * consts->ALPHA * consts->xi * consts->Fprime );
    gradz->y += U.x * (  .5 * consts->ALPHA * consts->Eprime );
    gradz->z += U.x * ( -.5 * ( 1. - consts->ALPHA ) * ( consts->sin / consts->R - consts->q * consts->Y11 * consts->cos ) - .5 * consts->ALPHA * consts->q * consts->Fprime );
  }
  if ( U.y ) {
    gradz->x += U.y * (  .5 * consts->ALPHA * consts->Eprime );
    gradz->y += U.y * (  .5 * ( 1. - consts->ALPHA ) * consts->ytilde * consts->X11 + .5 * consts->xi * consts->Y11 * consts->cos + .5 * consts->ALPHA * consts->eta * consts->Gprime );
    gradz->z += U.y * ( -.5 * ( 1. - consts->ALPHA ) * consts->dtilde * consts->X11 - .5 * consts->ALPHA * consts->q * consts->Gprime );
  }
  if ( U.z ) {
    gradz->x += U.z * (  .5 * ( 1. - consts->ALPHA ) * ( consts->sin / consts->R - consts->q * consts->Y11 * consts->cos ) - .5 * consts->ALPHA * consts->q * consts->Fprime );
    gradz->y += U.z * (  .5 * ( 1. - consts->ALPHA ) * consts->dtilde * consts->X11 - .5 * consts->ALPHA * consts->q * consts->Gprime );
    gradz->z += U.z * (  .5 * ( 1. - consts->ALPHA ) * ( consts->ytilde * consts->X11 + consts->xi * consts->Y11 * consts->cos ) + .5 * consts->ALPHA * consts->q * consts->Hprime );
  }
}
__inline void okada_add_grad_B( Matrix *grad, const OkadaConsts *consts, const Vector U ) {
  if ( U.x ) {
    grad->x.x += U.x * (  consts->xi * consts->xi * consts->q * consts->Y32 - ( 1. - consts->ALPHA ) / consts->ALPHA * consts->J1 * consts->sin );
    grad->x.y += U.x * (  consts->xi * consts->q / consts->R3 - ( 1. - consts->ALPHA ) / consts->ALPHA * consts->J2 * consts->sin );
    grad->x.z += U.x * ( -consts->xi * consts->q * consts->q * consts->Y32 - ( 1. - consts->ALPHA ) / consts->ALPHA * consts->J3 * consts->sin );
    grad->y.x += U.x * ( -consts->xi * consts->F - consts->dtilde * consts->X11 + ( 1. - consts->ALPHA ) / consts->ALPHA * ( consts->xi * consts->Y11 + consts->J4 ) * consts->sin );
    grad->y.y += U.x * ( -consts->E + ( 1. - consts->ALPHA ) / consts->ALPHA * ( 1. / consts->R + consts->J5 ) * consts->sin );
    grad->y.z += U.x * (  consts->q * consts->F - ( 1. - consts->ALPHA ) / consts->ALPHA * ( consts->q * consts->Y11 - consts->J6 ) * consts->sin );
    grad->z.x += U.x * ( -consts->xi * consts->Fprime - consts->ytilde * consts->X11 + ( 1. - consts->ALPHA ) / consts->ALPHA * consts->K1 * consts->sin );
    grad->z.y += U.x * ( -consts->Eprime + ( 1. - consts->ALPHA ) / consts->ALPHA * consts->ytilde * consts->D11 * consts->sin );
    grad->z.z += U.x * (  consts->q * consts->Fprime + ( 1. - consts->ALPHA ) / consts->ALPHA * consts->K2 * consts->sin );
  }
  if ( U.y ) {
    grad->x.x += U.y * ( consts->xi * consts->q / consts->R3 + ( 1. - consts->ALPHA ) / consts->ALPHA * consts->J4 * consts->sin * consts->cos );
    grad->x.y += U.y * ( consts->eta * consts->q / consts->R3 + consts->q * consts->Y11 + ( 1. - consts->ALPHA ) / consts->ALPHA * consts->J5 * consts->sin * consts->cos );
    grad->x.z += U.y * ( -consts->q * consts->q / consts->R3 + ( 1. - consts->ALPHA ) / consts->ALPHA * consts->J6 * consts->sin * consts->cos );
    grad->y.x += U.y * ( -consts->E + ( 1. - consts->ALPHA ) / consts->ALPHA * consts->J1 * consts->sin * consts->cos );
    grad->y.y += U.y * ( -consts->eta * consts->G - consts->xi * consts->Y11 * consts->sin + ( 1. - consts->ALPHA ) / consts->ALPHA * consts->J2 * consts->sin * consts->cos );
    grad->y.z += U.y * ( consts->q * consts->G + ( 1. - consts->ALPHA ) / consts->ALPHA * consts->J3 * consts->sin * consts->cos );
    grad->z.x += U.y * ( -consts->Eprime - ( 1. - consts->ALPHA ) / consts->ALPHA * consts->K3 * consts->sin * consts->cos );
    grad->z.y += U.y * ( -consts->eta * consts->Gprime - consts->xi * consts->Y11 * consts->cos - ( 1. - consts->ALPHA ) / consts->ALPHA * consts->xi * consts->D11 * consts->sin * consts->cos );
    grad->z.z += U.y * ( consts->q * consts->Gprime - ( 1. - consts->ALPHA ) / consts->ALPHA * consts->K4 * consts->sin * consts->cos );
  }
  if ( U.z ) {
    grad->x.x += U.z * ( -consts->xi * consts->q * consts->q * consts->Y32 - ( 1. - consts->ALPHA ) / consts->ALPHA * consts->J4 * consts->sin * consts->sin );
    grad->x.y += U.z * ( -consts->q * consts->q / consts->R3 - ( 1. - consts->ALPHA ) / consts->ALPHA * consts->J5 * consts->sin * consts->sin );
    grad->x.z += U.z * ( consts->q * consts->q * consts->q * consts->Y32 - ( 1. - consts->ALPHA ) / consts->ALPHA * consts->J6 * consts->sin * consts->sin );
    grad->y.x += U.z * ( consts->q * consts->F - ( 1. - consts->ALPHA ) / consts->ALPHA * consts->J1 * consts->sin * consts->sin );
    grad->y.y += U.z * ( consts->q * consts->G - ( 1. - consts->ALPHA ) / consts->ALPHA * consts->J2 * consts->sin * consts->sin );
    grad->y.z += U.z * ( -consts->q * consts->H - ( 1. - consts->ALPHA ) / consts->ALPHA * consts->J3 * consts->sin * consts->sin );
    grad->z.x += U.z * ( consts->q * consts->Fprime + ( 1. - consts->ALPHA ) / consts->ALPHA * consts->K3 * consts->sin * consts->sin );
    grad->z.y += U.z * ( consts->q * consts->Gprime + ( 1. - consts->ALPHA ) / consts->ALPHA * consts->xi * consts->D11 * consts->sin * consts->sin );
    grad->z.z += U.z * ( -consts->q * consts->Hprime + ( 1. - consts->ALPHA ) / consts->ALPHA * consts->K4 * consts->sin * consts->sin );
  }
}
__inline void okada_add_grad_C( Matrix *grad, const OkadaConsts *consts, const Vector U ) {
  if ( U.x ) {
    grad->x.x += U.x * (  ( 1. - consts->ALPHA ) * consts->Y0 * consts->cos - consts->ALPHA * consts->q * consts->Z0 );
    grad->x.y += U.x * ( -( 1. - consts->ALPHA ) * consts->xi * ( consts->cos / consts->R3 + 2. * consts->q * consts->Y32 * consts->sin ) + consts->ALPHA * 3. * consts->ctilde * consts->xi * consts->q / consts->R5 );
    grad->x.z += U.x * ( -( 1. - consts->ALPHA ) * consts->xi * consts->q * consts->Y32 * consts->cos + consts->ALPHA * consts->xi * ( 3. * consts->ctilde * consts->eta / consts->R5 - consts->z * consts->Y32 - consts->Z32 - consts->Z0 ) );
    grad->y.x += U.x * ( -( 1. - consts->ALPHA ) * consts->xi * consts->P * consts->cos - consts->ALPHA * consts->xi * consts->Q );
    grad->y.y += U.x * ( 2. * ( 1. - consts->ALPHA ) * ( consts->dtilde / consts->R3 - consts->Y0 * consts->sin ) * consts->sin - consts->ytilde / consts->R3 * consts->cos - consts->ALPHA * ( ( consts->ctilde + consts->dtilde ) / consts->R3 * consts->sin - consts->eta / consts->R3 - 3. * consts->ctilde * consts->ytilde * consts->q  / consts->R5 ) );
    grad->y.z += U.x * ( -( 1. - consts->ALPHA ) * consts->q / consts->R3 + ( consts->ytilde / consts->R3 - consts->Y0 * consts->cos ) * consts->sin + consts->ALPHA * ( ( consts->ctilde + consts->dtilde ) / consts->R3 * consts->cos + 3. * consts->ctilde * consts->dtilde * consts->q / consts->R5 - ( consts->Y0 * consts->cos + consts->q * consts->Z0 ) * consts->sin ) );
    grad->z.x += U.x * ( ( 1. - consts->ALPHA ) * consts->xi * consts->Pprime * consts->cos - consts->ALPHA * consts->xi * consts->Qprime );
    grad->z.y += U.x * ( 2. * ( 1. - consts->ALPHA ) * ( consts->ytilde / consts->R3 - consts->Y0 * consts->cos ) * consts->sin + consts->dtilde / consts->R3 * consts->cos - consts->ALPHA * ( ( consts->ctilde + consts->dtilde ) / consts->R3 * consts->cos + 3. * consts->ctilde * consts->dtilde * consts->q / consts->R5 ) );
    grad->z.z += U.x * ( ( consts->ytilde / consts->R3 - consts->Y0 * consts->cos ) * consts->cos - consts->ALPHA * ( ( consts->ctilde + consts->dtilde ) / consts->R3 * consts->sin - 3. * consts->ctilde * consts->ytilde * consts->q / consts->R5 - consts->Y0 * consts->sin * consts->sin + consts->q * consts->Z0 * consts->cos ) );
  }
  if ( U.y ) {
    grad->x.x += U.y * ( -( 1. - consts->ALPHA ) * consts->xi / consts->R3 * consts->cos + consts->xi * consts->q * consts->Y32 * consts->sin + consts->ALPHA * 3. * consts->ctilde * consts->xi * consts->q / consts->R5 );
    grad->x.y += U.y * ( -( 1. - consts->ALPHA ) * consts->ytilde / consts->R3 + consts->ALPHA * 3. * consts->ctilde * consts->eta * consts->q / consts->R5 );
    grad->x.z += U.y * ( consts->dtilde / consts->R3 - consts->Y0 * consts->sin + consts->ALPHA * consts->ctilde / consts->R3 * ( 1. - 3. * consts->q * consts->q / consts->R2 ) );
    grad->y.x += U.y * ( -( 1. - consts->ALPHA ) * consts->eta / consts->R3 + consts->Y0 * consts->sin * consts->sin - consts->ALPHA * ( ( consts->ctilde + consts->dtilde ) / consts->R3 * consts->sin - 3. * consts->ctilde * consts->ytilde * consts->q / consts->R5 ) );
    grad->y.y += U.y * ( ( 1. - consts->ALPHA ) * ( consts->X11 - consts->ytilde * consts->ytilde * consts->X32 ) - consts->ALPHA * consts->ctilde * ( ( consts->dtilde + 2. * consts->q * consts->cos ) * consts->X32 - consts->ytilde * consts->eta * consts->q * consts->X53 ) );
    grad->y.z += U.y * ( consts->xi * consts->P * consts->sin + consts->ytilde * consts->dtilde * consts->X32 + consts->ALPHA * consts->ctilde * ( ( consts->ytilde + 2. * consts->q * consts->sin ) * consts->X32 - consts->ytilde * consts->q * consts->q * consts->X53 ) );
    grad->z.x += U.y * ( -consts->q / consts->R3 + consts->Y0 * consts->sin * consts->cos - consts->ALPHA * ( ( consts->ctilde + consts->dtilde ) / consts->R3 * consts->cos + 3. * consts->ctilde * consts->dtilde * consts->q / consts->R5 ) );
    grad->z.y += U.y * ( ( 1. - consts->ALPHA ) * consts->ytilde * consts->dtilde * consts->X32 - consts->ALPHA * consts->ctilde * ( ( consts->ytilde - 2. * consts->q * consts->sin ) * consts->X32 + consts->dtilde * consts->eta * consts->q * consts->X53 ) );
    grad->z.z += U.y * ( -consts->xi * consts->Pprime * consts->sin + consts->X11 - consts->dtilde * consts->dtilde * consts->X32 - consts->ALPHA * consts->ctilde * ( ( consts->dtilde - 2. * consts->q * consts->cos ) * consts->X32 - consts->dtilde * consts->q * consts->q * consts->X53 ) );
  }
  if ( U.z ) {
    grad->x.x += U.z * ( ( 1. - consts->ALPHA ) * consts->xi / consts->R3 * consts->sin + consts->xi * consts->q * consts->Y32 * consts->cos + consts->ALPHA * consts->xi * ( 3. * consts->ctilde * consts->eta / consts->R5 - 2. * consts->Z32 - consts->Z0 ) );
    grad->x.y += U.z * ( ( 1. - consts->ALPHA ) * 2. * consts->Y0 * consts->sin - consts->dtilde / consts->R3 + consts->ALPHA * consts->ctilde / consts->R3 * ( 1. - 3. * consts->q * consts->q / consts->R2 ) );
    grad->x.z += U.z * ( -( 1. - consts->ALPHA ) * ( consts->ytilde / consts->R3 - consts->Y0 * consts->cos ) - consts->ALPHA * ( 3. * consts->ctilde * consts->eta * consts->q / consts->R5 - consts->q * consts->Z0 ) );
    grad->y.x += U.z * ( ( 1. - consts->ALPHA ) * ( consts->q / consts->R3 + consts->Y0 * consts->sin * consts->cos ) + consts->ALPHA * ( consts->z / consts->R3 * consts->cos + 3. * consts->ctilde * consts->dtilde * consts->q / consts->R5 - consts->q * consts->Z0 * consts->sin ) );
    grad->y.y += U.z * ( -( 1. - consts->ALPHA ) * 2. * consts->xi * consts->P * consts->sin - consts->ytilde * consts->dtilde * consts->X32 + consts->ALPHA * consts->ctilde * ( ( consts->ytilde + 2. * consts->q * consts->sin ) * consts->X32 - consts->ytilde * consts->q * consts->q * consts->X53 ) );
    grad->y.z += U.z * ( -( 1. - consts->ALPHA ) * ( consts->xi * consts->P * consts->cos - consts->X11 + consts->ytilde * consts->ytilde * consts->X32 ) + consts->ALPHA * consts->ctilde * ( ( consts->dtilde + 2. * consts->q * consts->cos ) * consts->X32 - consts->ytilde * consts->eta * consts->q * consts->X53 ) + consts->ALPHA * consts->xi * consts->Q );
    grad->z.x += U.z * ( -consts->eta / consts->R3 + consts->Y0 * consts->cos * consts->cos - consts->ALPHA * ( consts->z / consts->R3 * consts->sin - 3. * consts->ctilde * consts->ytilde * consts->q / consts->R5 - consts->Y0 * consts->sin * consts->sin + consts->q * consts->Z0 * consts->cos ) );
    grad->z.y += U.z * ( ( 1. - consts->ALPHA ) * 2. * consts->xi * consts->Pprime * consts->sin - consts->X11 + consts->dtilde * consts->dtilde * consts->X32 - consts->ALPHA * consts->ctilde * ( ( consts->dtilde - 2. * consts->q * consts->cos ) * consts->X32 - consts->dtilde * consts->q * consts->q * consts->X53 ) );
    grad->z.z += U.z * ( ( 1. - consts->ALPHA ) * ( consts->xi * consts->Pprime * consts->cos + consts->ytilde * consts->dtilde * consts->X32 ) + consts->ALPHA * consts->ctilde * ( ( consts->ytilde - 2. * consts->q * consts->sin ) * consts->X32 + consts->dtilde * consts->eta * consts->q * consts->X53 ) + consts->ALPHA * consts->xi * consts->Qprime );
  }
}

// +-----------------------+
// |                       |
// |       MAIN API        |
// |                       |
// +-----------------------+

typedef struct { // SourceParams
  double strike;
  double dip;
  double length;
  double width;
  double xbottom, ybottom, zbottom;
  double strikeslip, dipslip, opening;
} SourceParams;
Vector get_displacement( SourceParams *params, Vector *where, double poisson ) {
  const Angle dip = from_deg( params->dip );
  const Angle perp_strike = from_deg( params->strike - 90. );
  Vector U = { params->strikeslip / ( 2. * M_PI ), params->dipslip / ( 2. * M_PI ), params->opening / ( 2. * M_PI ) };
  Vector point = { where->x - params->xbottom, where->y - params->ybottom, where->z };
  OkadaConsts consts;
  int i;
  Vector disp = { 0 };
  if ( point.z == 0 ) {
    Vector_irot_xy( &point, perp_strike );
    for ( i = 0; i < 4; i++ ) {
      okada_get_consts( &consts, point, -params->zbottom, dip, params->length, params->width, i, poisson );
      if ( i & 1 ) Vector_flip( &U );
      okada_add_disp_B( &disp, &consts, U );
    }
    Vector_irot_yz( &disp, dip );
  }
  else {
    Vector v = { 0 };
    Vector_irot_xy( &point, perp_strike );
    for ( i = 0; i < 4; i++ ) {
      okada_get_consts( &consts, point, -params->zbottom, dip, params->length, params->width, i, poisson );
      if ( i & 1 ) Vector_flip( &U );
      okada_add_disp_A( &disp, &consts, U );
      okada_add_disp_B( &disp, &consts, U );
      okada_add_disp_C( &v, &consts, Vector_mul( U, point.z ) );
    }
    v.z = -v.z;
    point.z = -point.z;
    for ( i = 0; i < 4; i++ ) {
      okada_get_consts( &consts, point, -params->zbottom, dip, params->length, params->width, i, poisson );
      if ( i & 1 ) Vector_flip( &U );
      okada_add_disp_A( &disp, &consts, Vector_mul( U, -1.0 ) );
    }
    Vector_irot_yz( &disp, dip );
    Vector_irot_zy( &v, dip );
    Vector_iadd( &disp, v );
  }
  Vector_irot_yx( &disp, perp_strike );
  return disp;
}
Matrix get_gradient( SourceParams *params, Vector *where, double poisson ) {
  const Angle dip = from_deg( params->dip );
  const Angle perp_strike = from_deg( params->strike - 90. );
  Vector U = { params->strikeslip / ( 2. * M_PI ), params->dipslip / ( 2. * M_PI ), params->opening / ( 2. * M_PI ) };
  Vector point = { where->x - params->xbottom, where->y - params->ybottom, where->z };
  Matrix grad = { 0 };
  OkadaConsts consts;
  int i;
  if ( point.z == 0 ) {
    Vector h_z = { 0 };
    Vector_irot_xy( &point, perp_strike );
    for ( i = 0; i < 4; i++ ) {
      okada_get_more_consts( &consts, point, -params->zbottom, dip, params->length, params->width, i, poisson );
      if ( i & 1 ) Vector_flip( &U );
      okada_add_grad_Az( &(grad.z), &consts, Vector_mul( U,  2.0 ) );
      okada_add_grad_B( &grad, &consts, U );
      okada_add_disp_C( &h_z, &consts, U );
    }
    h_z.z = -h_z.z;
    Vector_irot_zy( &(h_z), dip );
    for ( i = 0; i < 3; i++ ) {
      Vector_irot_yz( (Vector *)&grad + i, dip );
    }
    Vector_iadd( &(grad.z), h_z );
  }
  else {
    Matrix h = {{ 0 }};
    Vector_irot_xy( &point, perp_strike );
    for ( i = 0; i < 4; i++ ) {
      okada_get_more_consts( &consts, point, -params->zbottom, dip, params->length, params->width, i, poisson );
      if ( i & 1 ) Vector_flip( &U );
      okada_add_grad_Axy( &grad, &consts, U );
      okada_add_grad_Az( &(grad.z), &consts, U );
      okada_add_grad_B( &grad, &consts, U );
      okada_add_grad_C( &h, &consts, Vector_mul( U, point.z ) );
      okada_add_disp_C( &(h.z), &consts, U );
    }
    point.z = -point.z;
    for ( i = 0; i < 4; i++ ) {
      okada_get_more_consts( &consts, point, -params->zbottom, dip, params->length, params->width, i, poisson );
      if ( i & 1 ) Vector_flip( &U );
      okada_add_grad_Axy( &grad, &consts, Vector_mul( U, -1.0 ) );
      okada_add_grad_Az( &(grad.z), &consts, U );
    }
    for ( i = 0; i < 3; i++ ) {
      ((Vector *)&h)[i].z = -((Vector *)&h)[i].z;
      Vector_irot_zy( (Vector *)&h + i, dip );
      Vector_irot_yz( (Vector *)&grad + i, dip );
    }
    Matrix_iadd( &grad, h );
  }
  Matrix_irot_xy( &grad, perp_strike );
  return grad;
}
EXPORT void get_displacements( Vector *out, SourceParams *params, Vector *where, double poisson, int count ) {
  while ( count-- ) Vector_copy( out++, get_displacement(params,where++,poisson) );
}
EXPORT void get_gradients( Matrix *out, SourceParams *params, Vector *where, double poisson, int count ) {
  while ( count-- ) Matrix_copy( out++, get_gradient(params,where++,poisson) );
}

// vim:foldmethod=syntax
