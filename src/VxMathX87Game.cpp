// BallanceMMO x87-soft headless engine (VX_X87_GAME): the 2002 VxMath.dll routines that run as x87 code, step for step
// on the IVP x87 model at the game thread's precision (IVP_X87_Control::current(): 24 bits before the process's first
// hull, 64 after) with nearest rounding.  Matrices are 16 floats, row-major [r][c].  Each one was compared bit for bit
// with the 2002 DLL called at 24 and at 64 bits of precision (2e5 / 1e5 inputs, special angles and singular matrices
// included).
#if defined(VX_X87_GAME)
#include <ivu_x87.hxx>

typedef IVP_U_Float80 F;

// VxMath.dll 0x1B280 Vx3DMatrixFromEulerAngles(VxMatrix&, float eax, float eay, float eaz):
// fcos/fsin of each angle stored to float, except cos(eaz) which stays on the stack (80-bit) until it is consumed by
// the m[0][0] product; |angle| <= 1e-8 (double constant) gives cos 1, sin 0 (cos z then the float constant 1.0).
void vx_x87_matrix_from_euler(float *m, float eax, float eay, float eaz)
{
    IVP_X87_Arithmetic x(IVP_X87_Control::current());
    const F eps = x.load(9.99999993922529e-09);   // the qword at 0x242a0238
    float cx = x.store_float(x.cosine(x.load(eax)));
    float cy = x.store_float(x.cosine(x.load(eay)));
    F cz = x.cosine(x.load(eaz));
    float sx = x.store_float(x.sine(x.load(eax)));
    float sy = x.store_float(x.sine(x.load(eay)));
    float sz = x.store_float(x.sine(x.load(eaz)));
    auto small = [&](float a) {
        const IVP_X87_Arithmetic::Relation r = x.compare(x.absolute(x.load(a)), eps);
        return r == IVP_X87_Arithmetic::LESS || r == IVP_X87_Arithmetic::EQUAL;   // test ah,0x41; jp
    };
    if (small(eax)) { cx = 1.0f; sx = 0.0f; }
    if (small(eay)) { cy = 1.0f; sy = 0.0f; }
    if (small(eaz)) { cz = x.load(1.0f); sz = 0.0f; }
    const float t24 = x.store_float(x.multiply(x.load(cx), cz));
    const float t28 = x.store_float(x.multiply(x.load(sz), x.load(cx)));
    const float t14 = x.store_float(x.multiply(x.load(sx), cz));
    const float t20 = x.store_float(x.multiply(x.load(sz), x.load(sx)));
    m[0] = x.store_float(x.multiply(cz, x.load(cy)));
    m[4] = x.store_float(x.subtract(x.multiply(x.load(t14), x.load(sy)), x.load(t28)));
    m[8] = x.store_float(x.add(x.multiply(x.load(t24), x.load(sy)), x.load(t20)));
    m[1] = x.store_float(x.multiply(x.load(sz), x.load(cy)));
    m[5] = x.store_float(x.add(x.multiply(x.load(t20), x.load(sy)), x.load(t24)));
    m[9] = x.store_float(x.subtract(x.multiply(x.load(t28), x.load(sy)), x.load(t14)));
    m[2] = x.store_float(x.negate(x.load(sy)));
    m[6] = x.store_float(x.multiply(x.load(sx), x.load(cy)));
    m[10] = x.store_float(x.multiply(x.load(cy), x.load(cx)));
    m[12] = m[13] = m[14] = 0.0f;
    m[3] = m[7] = m[11] = 0.0f;
    m[15] = 1.0f;
}

// VxMath.dll 0x1AD50 Vx3DInverseMatrix(VxMatrix &dest, const VxMatrix &m): the 3x3 cofactor inverse with the
// determinant's positive and negative terms summed apart on the stack (the first term through a float store), a
// zero determinant or |det/(pos-neg)| < 1e-15 leaving dest unwritten; the translation row from the stored result.
void vx_x87_inverse_matrix(float *d, const float *m)
{
    IVP_X87_Arithmetic x(IVP_X87_Control::current());
    auto L = [&](float v) { return x.load(v); };
    const F zero = L(0.0f);
    auto negative = [&](const F &v) { return x.compare(v, zero) == IVP_X87_Arithmetic::LESS; };
    const F p0 = x.multiply(L(m[5]), L(m[10]));
    const float t0 = x.store_float(p0);
    const float s = x.store_float(x.multiply(p0, L(m[0])));
    F pos, neg;
    if (negative(x.multiply(p0, L(m[0])))) { pos = L(0.0f); neg = L(s); }
    else { pos = L(s); neg = zero; }
    auto acc = [&](const F &t) { if (negative(t)) neg = x.add(neg, t); else pos = x.add(pos, t); };
    acc(x.multiply(x.multiply(L(m[6]), L(m[8])), L(m[1])));
    acc(x.multiply(x.multiply(L(m[2]), L(m[4])), L(m[9])));
    acc(x.negate(x.multiply(x.multiply(L(m[2]), L(m[8])), L(m[5]))));
    acc(x.negate(x.multiply(x.multiply(L(m[4]), L(m[10])), L(m[1]))));
    const F p4 = x.multiply(L(m[6]), L(m[9]));
    const float t4 = x.store_float(p4);
    acc(x.negate(x.multiply(p4, L(m[0]))));
    const F det = x.add(neg, pos);
    const float sdet = x.store_float(det);
    const IVP_X87_Arithmetic::Relation r = x.compare(det, zero);
    if (r == IVP_X87_Arithmetic::EQUAL) return;
    const F ratio = x.absolute(x.divide(L(sdet), x.subtract(pos, neg)));
    if (x.compare(ratio, x.load(1e-15)) == IVP_X87_Arithmetic::LESS) return;
    const F inv = x.divide(L(1.0f), L(sdet));
    auto cof = [&](int a, int b, int c, int e) { return x.subtract(x.multiply(L(m[a]), L(m[b])), x.multiply(L(m[c]), L(m[e]))); };
    d[0] = x.store_float(x.multiply(x.subtract(L(t0), L(t4)), inv));
    d[4] = x.store_float(x.negate(x.multiply(cof(4, 10, 6, 8), inv)));
    d[8] = x.store_float(x.multiply(cof(4, 9, 8, 5), inv));
    d[1] = x.store_float(x.negate(x.multiply(cof(10, 1, 2, 9), inv)));
    d[5] = x.store_float(x.multiply(cof(0, 10, 2, 8), inv));
    d[9] = x.store_float(x.negate(x.multiply(cof(0, 9, 8, 1), inv)));
    d[2] = x.store_float(x.multiply(cof(6, 1, 2, 5), inv));
    d[6] = x.store_float(x.negate(x.multiply(cof(6, 0, 2, 4), inv)));
    d[10] = x.store_float(x.multiply(cof(5, 0, 4, 1), inv));
    d[3] = d[7] = d[11] = 0.0f;
    d[12] = x.store_float(x.negate(x.add(x.add(x.multiply(L(m[13]), L(d[4])), x.multiply(L(d[0]), L(m[12]))),
                                         x.multiply(L(d[8]), L(m[14])))));
    d[13] = x.store_float(x.negate(x.add(x.add(x.multiply(L(d[5]), L(m[13])), x.multiply(L(d[1]), L(m[12]))),
                                         x.multiply(L(m[14]), L(d[9])))));
    d[14] = x.store_float(x.negate(x.add(x.add(x.multiply(L(d[10]), L(m[14])), x.multiply(L(m[13]), L(d[6]))),
                                         x.multiply(L(m[12]), L(d[2])))));
    d[15] = 1.0f;
}
#endif
