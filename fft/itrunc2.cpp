#define PD_RADIX_4_FORWARD_2X_ITRUNC2_OTRUNC4(X0, X1, X2, X3, Z0, Z1, Z2, Z3) \
{ \
    packed<double, 4> _x0, _x1, _x2, _x3, _y0, _y1, _y2, _y3; \
    packed<double, 4> _z0, _z1, _z2, _z3, _w0, _w1, _w2, _w3; \
    _x0.load(X0); \
        _z0.load(Z0); \
    _x0 = reduce_to_pm1n(_x0, n, ninv); \
        _z0 = reduce_to_pm1n(_z0, n, ninv); \
    _x1.load(X1); \
        _z1.load(Z1); \
    _y0 = _x0; \
        _w0 = _z0; \
    _y1 = _x1; \
        _w1 = _z1; \
    _y2 = _x0; \
        _w2 = _z0; \
    _y3 = _x1; \
        _w3 = _z1; \
    _y1 = mulmod2(_y1, w, n, ninv); \
        _w1 = mulmod2(_w1, w, n, ninv); \
    _y3 = mulmod2(_y3, iw, n, ninv); \
        _w3 = mulmod2(_w3, iw, n, ninv); \
    _x0 = add(_y0, _y1); \
        _z0 = add(_w0, _w1); \
    _x1 = sub(_y0, _y1); \
        _z1 = sub(_w0, _w1); \
    _x2 = add(_y2, _y3); \
        _z2 = add(_w2, _w3); \
    _x3 = sub(_y2, _y3); \
        _z3 = sub(_w2, _w3); \
    _x0.store(X0); \
        _z0.store(Z0); \
    _x1.store(X1); \
        _z1.store(Z1); \
    _x2.store(X2); \
        _z2.store(Z2); \
    _x3.store(X3); \
        _z3.store(Z3); \
}

#define PD_RADIX_4_FORWARD_2X_ITRUNC2_OTRUNC3(X0, X1, X2, X3, Z0, Z1, Z2, Z3) \
{ \
    packed<double, 4> _x0, _x1, _x2, _x3, _y0, _y1, _y2, _y3; \
    packed<double, 4> _z0, _z1, _z2, _z3, _w0, _w1, _w2, _w3; \
    _x0.load(X0); \
        _z0.load(Z0); \
    _x0 = reduce_to_pm1n(_x0, n, ninv); \
        _z0 = reduce_to_pm1n(_z0, n, ninv); \
    _x1.load(X1); \
        _z1.load(Z1); \
    _y0 = _x0; \
        _w0 = _z0; \
    _y1 = _x1; \
        _w1 = _z1; \
    _y2 = _x0; \
        _w2 = _z0; \
    _y3 = _x1; \
        _w3 = _z1; \
    _y1 = mulmod2(_y1, w, n, ninv); \
        _w1 = mulmod2(_w1, w, n, ninv); \
    _y3 = mulmod2(_y3, iw, n, ninv); \
        _w3 = mulmod2(_w3, iw, n, ninv); \
    _x0 = add(_y0, _y1); \
        _z0 = add(_w0, _w1); \
    _x1 = sub(_y0, _y1); \
        _z1 = sub(_w0, _w1); \
    _x2 = add(_y2, _y3); \
        _z2 = add(_w2, _w3); \
    _x0.store(X0); \
        _z0.store(Z0); \
    _x1.store(X1); \
        _z1.store(Z1); \
    _x2.store(X2); \
        _z2.store(Z2); \
}

#define PD_RADIX_4_FORWARD_2X_ITRUNC2_OTRUNC2(X0, X1, X2, X3, Z0, Z1, Z2, Z3) \
{ \
    packed<double, 4> _x0, _x1, _x2, _x3, _y0, _y1, _y2, _y3; \
    packed<double, 4> _z0, _z1, _z2, _z3, _w0, _w1, _w2, _w3; \
    _x0.load(X0); \
        _z0.load(Z0); \
    _x0 = reduce_to_pm1n(_x0, n, ninv); \
        _z0 = reduce_to_pm1n(_z0, n, ninv); \
    _x1.load(X1); \
        _z1.load(Z1); \
    _y0 = _x0; \
        _w0 = _z0; \
    _y1 = _x1; \
        _w1 = _z1; \
    _y2 = _x0; \
        _w2 = _z0; \
    _y3 = _x1; \
        _w3 = _z1; \
    _y1 = mulmod2(_y1, w, n, ninv); \
        _w1 = mulmod2(_w1, w, n, ninv); \
    _x0 = add(_y0, _y1); \
        _z0 = add(_w0, _w1); \
    _x1 = sub(_y0, _y1); \
        _z1 = sub(_w0, _w1); \
    _x0.store(X0); \
        _z0.store(Z0); \
    _x1.store(X1); \
        _z1.store(Z1); \
}


#define PD_RADIX_4_FORWARD_2X_ITRUNC2_OTRUNC1(X0, X1, X2, X3, Z0, Z1, Z2, Z3) \
{ \
    packed<double, 4> _x0, _x1, _x2, _x3, _y0, _y1, _y2, _y3; \
    packed<double, 4> _z0, _z1, _z2, _z3, _w0, _w1, _w2, _w3; \
    _x0.load(X0); \
        _z0.load(Z0); \
    _x0 = reduce_to_pm1n(_x0, n, ninv); \
        _z0 = reduce_to_pm1n(_z0, n, ninv); \
    _x1.load(X1); \
        _z1.load(Z1); \
    _y0 = _x0; \
        _w0 = _z0; \
    _y1 = _x1; \
        _w1 = _z1; \
    _y2 = _x0; \
        _w2 = _z0; \
    _y3 = _x1; \
        _w3 = _z1; \
    _y1 = mulmod2(_y1, w, n, ninv); \
        _w1 = mulmod2(_w1, w, n, ninv); \
    _x0 = add(_y0, _y1); \
        _z0 = add(_w0, _w1); \
    _x0.store(X0); \
        _z0.store(Z0); \
}



