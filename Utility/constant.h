// meson K, pi, phi, ks0: m0
const float kaonmass = 0.4937;
const float pionmass = 0.13957;
const float phimass = 1.0194;
const float ks0mass = 0.4976;
// p, lambda : m0
const float protonmass = 0.93827;
const float lambdamass = 1.11568;
const float casmass = 1.3217;
const float omgmass = 1.6724;

// number of constituent quark
const int nmeson = 2;
const int nbaryon = 3;

Float_t getNCQx(float pt, float m0, double n_mb)
{
    float mt_m0 = sqrt(pow(pt,2.0) + pow(m0,2.0)) - m0;
    float NCQx = mt_m0/n_mb;
    return NCQx;
}

Float_t getNCQy(float v2, double n_mb)
{
    float NCQy = v2/n_mb;
    return NCQy;
}
