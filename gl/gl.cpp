#include <gl/gl.hpp>
#include <sps/math.h>
#include <stddef.h>

#include <gl/gl_lut.hpp>
#include <gl/bessel_lut.hpp>

// TODO: Introduce dummy gl_lut.cpp file to support version without LUT

//
// Anonymous namespace for non-public functions
//
namespace {

  using gl::weights;
  using gl::abcissas;

  // TODO: Make inline, change Mathematica script to order weight and values equally
  gl::GLNode GLPairTabulated ( size_t l, size_t k )
  {
    gl::GLNode retval = gl::GLNode();
    size_t half = (l-1)/2;
    if ((l > 1) && (k < l)) {
      // Odd, verified
      if (l & 1) {
        if (k == half) {
          retval.value   = 0.0;
          retval.weight  = weights[l-2][0];
        } else if (k > half) {
          retval.value  = abcissas[l-2][2*half-k];
          // Weight are opposite
          retval.weight = weights[l-2][k-half];
        } else {
          retval.value  = -abcissas[l-2][k];
          // Weight are opposite
          retval.weight = weights[l-2][half-k];
        }
      }
      // Even, verified
      else {
        if (k > half) {
          // 2,3, half = 1
          retval.value  = abcissas[l-2][1+2*half-k];
          retval.weight = weights[l-2][k-half-1];
        } else {
          // 0,1
          retval.value  = -abcissas[l-2][k];
          retval.weight = weights[l-2][half-k];
        }
      }
    }
    return retval;
  }
}

namespace gl {
  double GLQuad(size_t n, double (*f)(double,void*), void* data, double a, double b)
  {
    const double m = 0.5 * (b - a);
    const double p = 0.5 * (a + b);
    double quad = 0.0;
    double x,w;
    for (size_t i = 0 ; i < n ; i++) {
      gl::GLNode qp = GL(n,i);
      x = qp.value;
      w = qp.weight;
      double s = m * x + p;
      quad += w*f(s,data);
    }
    return m*quad;
  }

  // Only function that needs to be public
  gl::GLNode GL(size_t l, size_t k)
  {

    // Inverted, 0.99, 0.5,.... -0.99
    if (l < _GL_LUT_TABLE_SIZE) {
      return GLPairTabulated(l,k);
    } else {
      if ( l < ( 2 * (k+1) - 1 ) ) {
        return GLS( l, l - k - 1 );
      } else {
        gl::GLNode P = GLS( l, k );
        P.value = - P.value;
        return P;
      }
    }
  }

  // TODO: Make inline
  gl::GLNode GLS(size_t n, size_t k)
  {

    //
    //  First get the Bessel zero
    //
    double w = 1.0/(n+0.5);
    double nu = besseljzero((int)k);
    double theta = w*nu;
    double x = theta*theta;
    //
    //  Get the asymptotic BesselJ(1,nu) squared
    //
    double B = besselj1squared((int)k);
    //
    //  Get the Chebyshev interpolants for the nodes... (6th order)
    //
    double SF1T = (((((-1.29052996274280508473467968379e-12*x +2.40724685864330121825976175184e-10)*x -3.13148654635992041468855740012e-8)*x +0.275573168962061235623801563453e-5)*x -0.148809523713909147898955880165e-3)*x +0.416666666665193394525296923981e-2)*x -0.416666666666662959639712457549e-1;
    double SF2T = (((((+2.20639421781871003734786884322e-9*x  -7.53036771373769326811030753538e-8)*x  +0.161969259453836261731700382098e-5)*x -0.253300326008232025914059965302e-4)*x +0.282116886057560434805998583817e-3)*x -0.209022248387852902722635654229e-2)*x +0.815972221772932265640401128517e-2;
    double SF3T = (((((-2.97058225375526229899781956673e-8*x  +5.55845330223796209655886325712e-7)*x  -0.567797841356833081642185432056e-5)*x +0.418498100329504574443885193835e-4)*x -0.251395293283965914823026348764e-3)*x +0.128654198542845137196151147483e-2)*x -0.416012165620204364833694266818e-2;
    //
    //  ...and for the weights (9th order)
    //
    double WSF1T = ((((((((-2.20902861044616638398573427475e-14*x+2.30365726860377376873232578871e-12)*x-1.75257700735423807659851042318e-10)*x+1.03756066927916795821098009353e-8)*x-4.63968647553221331251529631098e-7)*x+0.149644593625028648361395938176e-4)*x-0.326278659594412170300449074873e-3)*x+0.436507936507598105249726413120e-2)*x-0.305555555555553028279487898503e-1)*x+0.833333333333333302184063103900e-1;
    double WSF2T = (((((((+3.63117412152654783455929483029e-12*x+7.67643545069893130779501844323e-11)*x-7.12912857233642220650643150625e-9)*x+2.11483880685947151466370130277e-7)*x-0.381817918680045468483009307090e-5)*x+0.465969530694968391417927388162e-4)*x-0.407297185611335764191683161117e-3)*x+0.268959435694729660779984493795e-2)*x-0.111111111111214923138249347172e-1;
    double WSF3T = (((((((+2.01826791256703301806643264922e-9*x-4.38647122520206649251063212545e-8)*x+5.08898347288671653137451093208e-7)*x-0.397933316519135275712977531366e-5)*x+0.200559326396458326778521795392e-4)*x-0.422888059282921161626339411388e-4)*x-0.105646050254076140548678457002e-3)*x-0.947969308958577323145923317955e-4)*x+0.656966489926484797412985260842e-2;
    //
    //  Then refine with the paper expansions.
    //
    double NuoSin = nu/sin(theta);
    double BNuoSin = B*NuoSin;
    double WInvSinc = w*w*NuoSin;
    double WIS2 = WInvSinc*WInvSinc;
    //
    //  Finally compute the node and the weight.
    //
    theta = w*(nu + theta * WInvSinc * (SF1T + WIS2*(SF2T + WIS2*SF3T)));
    double Deno = BNuoSin + BNuoSin * WIS2*(WSF1T + WIS2*(WSF2T + WIS2*WSF3T));
    double weight = ( 2.0 * w ) / Deno;

    gl::GLNode P;
    P.value = cos(theta);
    P.weight = weight;

    return P;
  }

// TODO: Make inline
  double besseljzero(int k)
  {
    if(k > _BESSELJZERO_LUT_TABLE_SIZE - 1) {
      double z = M_PI*(k + 1.0 - 0.25);
      double r = 1.0/z;
      double r2 = r*r;
      z = z + r*(0.125+r2*(-0.807291666666666666666666666667e-1+r2*(0.246028645833333333333333333333+r2*(-1.82443876720610119047619047619+r2*(25.3364147973439050099206349206+r2*(-567.644412135183381139802038240+r2*(18690.4765282320653831636345064+r2*(-8.49353580299148769921876983660e5+5.09225462402226769498681286758e7*r2))))))));
      return z;
    } else {
      return JZ[k];
    }
  }

// TODO: Make inline
  double besselj1squared(int k)
  {
    if(k > _BESSELJ_1_SQUARED_LUT_TABLE_SIZE - 1) {
      double x = 1.0/(k + 1.0 - 0.25);
      double x2 = x*x;
      return x * (0.202642367284675542887758926420 + x2*x2*(-0.303380429711290253026202643516e-3 + x2*(0.198924364245969295201137972743e-3 + x2*(-0.228969902772111653038747229723e-3+x2*(0.433710719130746277915572905025e-3+x2*(-0.123632349727175414724737657367e-2+x2*(0.496101423268883102872271417616e-2+x2*(-0.266837393702323757700998557826e-1+.185395398206345628711318848386*x2))))))));
    } else {
      return J1[k];
    }
  }
}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
